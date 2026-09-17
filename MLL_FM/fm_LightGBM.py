import pandas as pd
import numpy as np
import shap
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import uproot
import lightgbm as lgb
import json
import ROOT

train_fraction = 0.8
num_boost_round = 2000
early_stopping_rounds = 50
num_leaves = 15
learning_rate = 0.05
min_data_in_leaf = 400
feature_fraction = 0.8
bagging_fraction = 0.8
bagging_freq = 1
lambda_l2 = 1.0

def add_relative_features(df):

    df = df.copy()
    g = df.groupby("event_id", group_keys=False)
    
    df["nll_rank"] = g["nll_weighted"].rank(
        ascending=True,
        method="first"
    )

    df["total_pe_rank"] = g["total_pe"].rank(
        ascending=False,
        method="first"
    )

    df["nhits_rank"] = g["nhits"].rank(
        ascending=False,
        method="first"
    )

    df["close_totalpe_rank"] = g["close_totalpe"].rank(
        ascending=False,
        method="first"
    )

    df["close_nhits_rank"] = g["close_nhits"].rank(
        ascending=False,
        method="first"
    )

    df["abs_dt_nearest_flash_rank"] = g["abs_dt_nearest_flash"].rank(
        ascending=True,
        method="first"
    )

    return df

if __name__ == "__main__":
    with open("./configs/ana_config.json", "r") as f:
        ana_config = json.load(f)
    sample_config_file = ana_config["sample_config_file"]
    use_preselection   = ana_config["use_preselection"]
    calib_method       = ana_config["calib_method"]

    with open("./configs/"+sample_config_file, "r") as f:
        sample_config = json.load(f)
    input_dir       = sample_config["input_dir"]
    geom_identifier = sample_config["geom_identifier"]

    with open("./configs/"+geom_identifier+".json", "r") as f:
        dune_geom_config = json.load(f)
    anode_x = dune_geom_config["anode_x"]
    print(f"Using geometry: {geom_identifier} with anode_x = {anode_x} cm")

    presel_suf = "_preselected" if use_preselection else ""
    df = uproot.open(input_dir+"MLL_Features_"+geom_identifier+"_"+calib_method+presel_suf+".root")["feature_tree"].arrays(library="pd")
    # print if any nan
    print("NaN values in the dataset:")
    print(df.isna().sum())
    
    unique_events = df["event_id"].unique()
    n_train       = int(len(unique_events) * train_fraction)
    train_events  = unique_events[:n_train]
    val_events    = unique_events[n_train:]

    train_df = df[df["event_id"].isin(train_events)].copy()
    val_df   = df[df["event_id"].isin(val_events)].copy()

    train_df = train_df.sort_values("event_id")
    val_df   = val_df.sort_values("event_id")

    print("Before adding relative features:")
    train_df = add_relative_features(train_df)
    print("After adding relative features")
    val_df   = add_relative_features(val_df)
    print("After adding relative features")

    # save val_df to a csv file for later use
    val_df.to_csv(input_dir+"val_df_"+geom_identifier+"_"+calib_method+presel_suf+".csv", index=False)

    print("Train events:", len(train_df["event_id"].unique()))
    print("Val events:", len(val_df["event_id"].unique()))
    print("Overlap:", set(train_events) & set(val_events))

    features = [
        "nll",
        "nll_weighted",
        "reco_term_mean", "reco_term_std", "reco_term_max", "reco_term_min",
        "noreco_term_mean", "noreco_term_std", "noreco_term_max", "noreco_term_min",
        "exp_ph_sum", "nhit_expected",
        "time_diff",
        "exp_reco_ratio",
        "nhits_expnhits_ratio",
        "exp_close_totalpe_ratio",
        "close_nhits_exp_ratio",
        "exp_near_totalpe_ratio",
        "charge", "max_charge",
        "y_reco", "z_reco",
        "charge_nhits",
        "total_pe", "max_pe",
        "totalpe_nhits_ratio",
        "nhits", "flash_reco_y", "flash_reco_z",
        "n_close_flashes",
        "close_totalpe",
        "close_nhits",
        "dt_nearest_flash",
        "abs_dt_nearest_flash",
        "n_near_flashes",
        "near_totalpe",
        "near_nhits",
        "e_reco",
        "flash_cluster_dist",
        "nll_rank",
        "total_pe_rank",
        "nhits_rank",
        "close_totalpe_rank",
        "close_nhits_rank",
        "abs_dt_nearest_flash_rank"
    ]


    train_data = lgb.Dataset(train_df[features], label=train_df["my_label"], group=train_df.groupby("event_id").size().to_list())
    val_data   = lgb.Dataset(val_df[features],   label=val_df["my_label"],   group=val_df.groupby("event_id").size().to_list())

    params = {
            "objective": "lambdarank",
            "metric": "ndcg",
            "ndcg_eval_at": [1, 3],
            "learning_rate": learning_rate,
            "num_leaves": num_leaves,
            "min_data_in_leaf": min_data_in_leaf,
            "feature_fraction": feature_fraction,
            "bagging_fraction": bagging_fraction,
            "bagging_freq": bagging_freq,
            "lambda_l2": lambda_l2,
            "verbose": -1
    }

    evals_result = {}

    model = lgb.train(
        params, 
        train_data,
        valid_sets=[train_data, val_data],
        valid_names=["train", "val"],
        num_boost_round=num_boost_round,
        callbacks=[
            lgb.early_stopping(early_stopping_rounds),
            lgb.record_evaluation(evals_result)
        ]
    )

    val_df["score"] = model.predict(val_df[features])
    best = val_df.loc[val_df.groupby("event_id")["score"].idxmax()]
    top1_acc = best["my_label"].mean()
    top1_acc_2 = len(best[best["my_label"] == 2]) / len(best)
    top1_acc_1_or_2 = len(best[best["my_label"] >= 1]) / len(best)
    print(f"Top-1 Accuracy: {top1_acc:.4f} (it makes sense only if bynary labels are used)")
    print(f"Top-1 Accuracy (label 2): {top1_acc_2:.4f}")
    print(f"Top-1 Accuracy (label 1 or 2): {top1_acc_1_or_2:.4f}")
    # Counts the fraction of mismatch relative to cheating: namely, the fraction of events where the best flash (according to the model)
    # has label 0, but there was at least a flash with label 1 or 2 in the same event.
    cheat_mismatch = 0
    model_mismatch = 0
    tryes = 0
    for event_id, group in val_df.groupby("event_id"):
        tryes += 1
        best_flash = group.loc[group["score"].idxmax()]
        if best_flash["my_label"] == 0:
            model_mismatch += 1
        if not (group["my_label"] >= 1).any():
            cheat_mismatch += 1
    print(f"Cheat Efficiency: {tryes - cheat_mismatch} / {tryes} = {(tryes - cheat_mismatch)/tryes:.4f}")
    print(f"Model Efficiency: {tryes - model_mismatch} / {tryes} = {(tryes - model_mismatch)/tryes:.4f}")
    print(f"Model/Cheat mismatches: {model_mismatch/cheat_mismatch:.4f}")

    
    # --- PLOTS ---------------------------------------------------------------
    # --- Prepare data ---
    print("Preparing data for plots...")
    train_ndcg = evals_result["train"]["ndcg@1"]
    val_ndcg   = evals_result["val"]["ndcg@1"]

    importance = model.feature_importance(importance_type="gain")
    feature_names = model.feature_name()

    # sort importance
    sorted_idx = np.argsort(importance)[::-1]
    top_k = min(10, len(importance))
    top_idx = sorted_idx[:top_k]

    # SHAP
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(val_df[features])

    # --- Create canvas ---
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # =========================================================
    # (1) Training curves
    # =========================================================
    print("Plotting training curves...")
    ax = axes[0, 0]
    ax.plot(train_ndcg, label="Train")
    ax.plot(val_ndcg, label="Validation")
    ax.set_xlabel("Iteration")
    ax.set_ylabel("NDCG@1")
    ax.set_title("Training Curve")
    ax.legend()

    # =========================================================
    # (2) Feature importance
    # =========================================================
    print("Plotting feature importance...")
    ax = axes[0, 1]
    vals = importance[top_idx][::-1]
    names = np.array(feature_names)[top_idx][::-1]

    ax.barh(range(len(vals)), vals)
    ax.set_yticks(range(len(vals)))
    ax.set_yticklabels(names)
    ax.set_yticks(range(top_k))
    ax.set_yticklabels(np.array(feature_names)[top_idx][::-1])
    ax.set_title("Feature Importance (Gain)")
    ax.set_xlabel("Gain")

    # =========================================================
    # (3) Score distribution
    # =========================================================
    print("Plotting score distribution...")
    ax = axes[1, 0]
    using_lable_2 = (val_df["my_label"] == 2).any()
    if using_lable_2:
        sns.histplot(
                val_df[val_df["my_label"] == 2]["score"],
                label="Target",
                kde=True,
                stat="density",
                ax=ax
        )
    sns.histplot(
            val_df[val_df["my_label"] == 1]["score"],
            label= "Purity > 0" if using_lable_2 else "True",
            kde=True,
            stat="density",
            ax=ax
    )
    sns.histplot(
            val_df[val_df["my_label"] == 0]["score"],
            label="False",
            kde=True,
            stat="density",
            ax=ax
    )
    ax.set_title("Score Distribution")
    ax.set_xlabel("Score")
    ax.legend()

    # =========================================================
    # (4) SHAP summary (bar only)
    # =========================================================
    print("Plotting SHAP summary...")
    ax = axes[1, 1]
    mmm = 7
    if mmm == 7: # to avoid Code is unreachable
        shap.summary_plot(
                shap_values,
                val_df[features],
                plot_type="bar",
                show=False
        )

    # SHAP ignores axes, so we force placement
    plt.sca(ax)

    # =========================================================
    # Final layout
    # =========================================================
    plt.tight_layout()
    plt.show()
    fig.savefig(input_dir+"MLL_LightGBM_plots_"+geom_identifier+"_"+calib_method+presel_suf+".png")

    # -------------------------------------------------------------------------
    # -------------------------------------------------------------------------

    # --- DEBUGS ---------------------------------------------------------------
    val_df["score"] = model.predict(val_df[features])
    best = val_df.loc[val_df.groupby("event_id")["score"].idxmax()]
    failures = best[best["my_label"] == 0]
    print("Failures:")
    print(failures.head())

    event_id = failures.iloc[0]["event_id"]
    group = val_df[val_df["event_id"] == event_id]
    print(f"Event ID: {event_id}")
    print(group.sort_values("score", ascending=False))

    # Create TEfficiency objects to store the mathcing efficiency and a function of x_true
    # of range (0, x_max) with 30 bins. Store it in a root file.
    out_file = ROOT.TFile(input_dir+"/MLL_LightGBM_efficiency_"+geom_identifier+"_"+calib_method+presel_suf+".root", "RECREATE")
    out_file.cd()
    max_drift = df["x_true"].max().max() if geom_identifier == "dune10kt" else df["x_true"].max().max() - df["x_true"].min().min()
    print(f"Max drift distance: {max_drift} cm")
    he_model  = ROOT.TEfficiency("he_eff_drift_model", "Matching Efficiency;Drift [cm];Efficiency", 30, 0., max_drift)
    he_max_pe = ROOT.TEfficiency("he_eff_drift_maxpe", "Matching Efficiency (max_pe);Drift [cm];Efficiency", 30, 0, max_drift)
    he_nll    = ROOT.TEfficiency("he_eff_drift_nll", "Matching Efficiency (nll);Drift [cm];Efficiency", 30, 0, max_drift)
    he_cheat = ROOT.TEfficiency("he_eff_drift_cheat", "Matching Efficiency (cheat);Drift [cm];Efficiency", 30, 0, max_drift)
    # Create a dictionary of TEfficiency for e_true in the range (5,7) (7,9) (9,11) (11,13) (13,15) (15,17)
    # givin the mid-point as name (6, 8, 10, 12, 14, 16)
    he_model_dict = {}
    he_max_pe_dict = {}
    he_nll_dict = {}
    he_cheat_dict = {}
    energy_bins = [(5,7), (7,9), (9,11), (11,13), (13,15), (15,17)]
    for e_true_min, e_true_max in energy_bins:
        e_true_mid = int((e_true_min + e_true_max) / 2)
        he_model_dict[e_true_mid] = ROOT.TEfficiency(f"he_eff_drift_model_e_{e_true_mid}", f"Matching Efficiency (e_true in ({e_true_min},{e_true_max}));Drift [cm];Efficiency", 30, 0., max_drift)
        he_max_pe_dict[e_true_mid] = ROOT.TEfficiency(f"he_eff_drift_maxpe_e_{e_true_mid}", f"Matching Efficiency (max_pe, e_true in ({e_true_min},{e_true_max}));Drift [cm];Efficiency", 30, 0, max_drift)
        he_nll_dict[e_true_mid] = ROOT.TEfficiency(f"he_eff_drift_nll_e_{e_true_mid}", f"Matching Efficiency (nll, e_true in ({e_true_min},{e_true_max}));Drift [cm];Efficiency", 30, 0, max_drift)
        he_cheat_dict[e_true_mid] = ROOT.TEfficiency(f"he_eff_drift_cheat_e_{e_true_mid}", f"Matching Efficiency (cheat, e_true in ({e_true_min},{e_true_max}));Drift [cm];Efficiency", 30, 0, max_drift)
    
    for _, row in best.iterrows():
        catch = 1 if int(row["my_label"] > 0) else 0
        x_drift = abs(row["x_true"]) if geom_identifier == "dune10kt" else anode_x - row["x_true"]
        he_model.Fill(catch, x_drift)
        e_true = row["e_true"]
        energy_bin = None
        for e_true_min, e_true_max in energy_bins:
            if e_true_min <= e_true < e_true_max:
                energy_bin = (e_true_min, e_true_max)
                break
        if energy_bin is not None:
            e_true_mid = (energy_bin[0] + energy_bin[1]) / 2
            he_model_dict[e_true_mid].Fill(catch, x_drift)
        
        # catch = 1 if row["purity"] > 0. else 0
        # he_model.Fill(catch, abs(row["x_true"]))



    # df.loc[:, "my_label"] = (df["purity"] > 0).astype(int)
    tryes2 = 0
    max_pe_efficiency = 0
    max_nll_weighted_efficiency = 0
    cheat_efficiency = 0
    for event_id, group in df.groupby("event_id"):
        best_max_pe = group.loc[group["total_pe"].idxmax()]
        best_max_nll_weighted = group.loc[group["nll_weighted"].idxmin()]
        best_cheat = group.loc[group["purity"].idxmax()]

        catch_max_pe = 1 if int(best_max_pe["my_label"]) > 0 else 0
        catch_max_nll_weighted = 1 if int(best_max_nll_weighted["my_label"]) > 0 else 0
        catch_cheat = 1 if best_cheat["purity"] > 0 else 0

        # catch_max_pe = 1 if best_max_pe["purity"] > 0. else 0
        # catch_max_nll_weighted = 1 if best_max_nll_weighted["purity"] > 0. else 0
        # catch_cheat = 1 if best_cheat["purity"] > 0. else 0

        x_drift = abs(best_max_pe["x_true"]) if geom_identifier == "dune10kt" else anode_x - best_max_pe["x_true"]
        he_max_pe.Fill(catch_max_pe, x_drift)
        he_nll.Fill(catch_max_nll_weighted, x_drift)
        he_cheat.Fill(catch_cheat, x_drift)

        tryes2 += 1
        if catch_max_pe == 1:
            max_pe_efficiency += 1
        if catch_max_nll_weighted == 1:
            max_nll_weighted_efficiency += 1
        if catch_cheat == 1:
            cheat_efficiency += 1

        energy_bin = None
        e_true = best_max_pe["e_true"]
        for e_true_min, e_true_max in energy_bins:
            if e_true_min <= e_true < e_true_max:
                energy_bin = (e_true_min, e_true_max)
                break
        if energy_bin is not None:
            e_true_mid = (energy_bin[0] + energy_bin[1]) / 2
            he_max_pe_dict[e_true_mid].Fill(catch_max_pe, x_drift)
            he_nll_dict[e_true_mid].Fill(catch_max_nll_weighted, x_drift)
            he_cheat_dict[e_true_mid].Fill(catch_cheat, x_drift)

    print(f"Max PE Efficiency: {max_pe_efficiency} / {tryes2} = {max_pe_efficiency/tryes2:.4f}")
    print(f"Max NLL Weighted Efficiency: {max_nll_weighted_efficiency} / {tryes2} = {max_nll_weighted_efficiency/tryes2:.4f}")
    print(f"Cheat Efficiency: {cheat_efficiency} / {tryes2} = {cheat_efficiency/tryes2:.4f}")

    # Save all the "performance print" into a text file
    with open(input_dir+"/MLL_LightGBM_performance_"+geom_identifier+"_"+calib_method+"_"+presel_suf+".txt", "w") as f:
        f.write(f"Top-1 Accuracy: {top1_acc:.4f} (it makes sense only if bynary labels are used)\n")
        f.write(f"Top-1 Accuracy (label 2): {top1_acc_2:.4f}\n")
        f.write(f"Top-1 Accuracy (label 1 or 2): {top1_acc_1_or_2:.4f}\n")
        f.write(f"Cheat Efficiency: {(tryes - cheat_mismatch)}/{tryes} = {(tryes - cheat_mismatch)/tryes:.4f}\n")
        f.write(f"Model Efficiency: {(tryes - model_mismatch)}/{tryes} = {(tryes - model_mismatch)/tryes:.4f}\n")
        f.write(f"Model/Cheat mismatches: {model_mismatch}/{cheat_mismatch} = {model_mismatch/cheat_mismatch:.4f}\n")
        f.write(f"Max PE Efficiency: {max_pe_efficiency}/{tryes2} = {max_pe_efficiency/tryes2:.4f}\n")
        f.write(f"Max NLL Weighted Efficiency: {max_nll_weighted_efficiency}/{tryes2} = {max_nll_weighted_efficiency/tryes2:.4f}\n")
        f.write(f"Cheat Efficiency: {cheat_efficiency}/{tryes2} = {cheat_efficiency/tryes2:.4f}\n")

    he_model.Write()
    he_max_pe.Write()
    he_nll.Write()
    he_cheat.Write()
    for e_true_mid in he_model_dict:
        he_model_dict[e_true_mid].Write()
        he_max_pe_dict[e_true_mid].Write()
        he_nll_dict[e_true_mid].Write()
        he_cheat_dict[e_true_mid].Write()
    out_file.Close()
