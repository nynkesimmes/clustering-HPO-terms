# -*- coding: utf-8 -*-
import os
from phenopy.build_hpo import generate_annotated_hpo_network
from phenopy.score import Scorer
import pandas as pd
import numpy as np
from scipy import stats
from imblearn.under_sampling import RandomUnderSampler

from hpoterms_dataframe import update_df
from similarity_measures import create_sem_sim_matrix
from similarity_measures import create_graph_sim_matrix_ged_gm
from similarity_measures import create_graph_sim_matrix_mcs
from clustering import cluster_spectral
from clustering import cluster_kmeans


"""
Created on Tue Nov 16 15:45:10 2021

@author: nynkesimmes
"""


def prepare_data():
    """
    Load all all necessary HPO files and call function that updates data.

    Returns
    -------
    df_updated : TYPE
        DESCRIPTION.
    hpo_network : TYPE
        DESCRIPTION.

    """
    # generate the hpo network and supporting objects

    # data directory; instantiate path to phenopy folder
    phenopy_data_directory = os.path.join(
        os.getenv('CONDA_PREFIX'), '.phenopy/data')
    # print(os.path.join(os.environ.get('CONDA_PREFIX')))

    # files used in building the annotated HPO network
    obo_file_path = os.path.join(phenopy_data_directory, 'hp.obo')
    disease_to_phenotype_file_path = os.path.join(
        phenopy_data_directory, 'phenotype.hpoa')

    hpo_network, alt2prim, disease_records = \
        generate_annotated_hpo_network(obo_file_path,
                                       disease_to_phenotype_file_path,)

    # load data and update data such that it corresponds with hpo.obo
    # make sure that the file is in this folder or use correct path
    # df_total = pd.read_pickle(r'C:\Users\nynke\Google Drive\BSc thesis\Data and code\BSc thesis 1\data_total.pickle')
    df_total = pd.read_pickle('data_total.pickle')
    df_updated = update_df(df_total, obo_file_path)  # update data
    return df_updated, hpo_network


def pipeline_single_run(df_total_updated, hpo_network, data_sample, sim_measures, measures_combi, cluster_method):
    """
    Apply whole method to total data set using helper functions.

    Parameters
    ----------
    df : DataFrame
        Updated total dataframe with 533 patients.
    hpo_network : classes.multidigraph.MultiDiGraph
        HPO network.
    data_sample : str
        String that indicates the concerned data set.
    sim_measures : tuple
        Strings of all the measures that will be applied.
    measures_combi : tuple
        Strings of all the combinations that will be applied.
    cluster_method : string
        String that indicates the clustering algorithm that will be applied.

    Returns
    -------
    results : list
        List of lists containing:
            - the string of the concerned data set,
            - the applied (combination of) measure(s)
            - and the resulting ARI

    """
    results = []
    for sim_measure in sim_measures:
        # # Create similarity matrices; takes some time.
        # # Once done, matrices can be loaded from file explorer
        # create_sim_matrices(df, hpo_network, data_sample, sim_measure, 0)

        # Apply clustering; nr_run = 0, because only 1 run for total data set
        ARI = cluster_sim_matrices(df_total_updated, data_sample, sim_measure, 0, cluster_method)
        results.append([data_sample, sim_measure, ARI])

    for combi in measures_combi:
        ARI = cluster_sim_matrices(df_total_updated, data_sample, combi, 0, cluster_method)
        results.append([data_sample, combi, ARI])

    return results


def pipeline_mtp_runs(df_updated_total, hpo_network, data_sample, sim_measures, measures_combi, cluster_method, total_nr_runs):
    """
    Apply whole method multiple times to undersampled data set using helper functions.

    Parameters
    ----------
    df : DataFrame
        Updated total dataframe with 533 patients.
    hpo_network : classes.multidigraph.MultiDiGraph
        HPO network.
    data_sample : str
        String that indicates the concerned data set.
    sim_measures : tuple
        Strings of all the measures that will be applied.
    measures_combi : tuple
        Strings of all the combinations that will be applied.
    cluster_method : string
        String that indicates the clustering algorithm that will be applied.
    total_nr_runs : int
        Integer that indicates number of runs that need to be done.

    Returns
    -------
    results : list
        List of lists containing:
            - the string of the concerned data set,
            - the applied (combination of) measure(s)
            - and the resulting ARI
        and for the individually applied measures also:
            - list of all obtained ARIs
            - stats

    """
    # # THIS PART BELOW (CREATE DATA SETS AND MATRICES) TAKES A LOT OF TIME.
    # # ONCE DONE, DATA SETS AND MATRICES CAN BE LOADED FROM FILE EXPLORER.
    # # SO THE PART BELOW CAN BE COMMENTED

    # rus = RandomUnderSampler()  # create instance of RandomUnderSampler

    # # Create and save multiple data sets
    # for nr_run in range(total_nr_runs):  # run multiple times
    #     if (data_sample == 'undersampled'):
    #         df, corresp_labels = rus.fit_resample(
    #             df_updated_total.iloc[:, :6],
    #             df_updated_total['label'])

    #     if (data_sample == 'undersampled_withoutspop'):
    #         df_updated_withoutspop = df_updated_total[7:553]  # create data set without SPOP_1 and SPOP_2
    #         df, corresp_labels = rus.fit_resample(
    #             df_updated_withoutspop.iloc[:, :6],
    #             df_updated_withoutspop['label'])

    #     # Save every created data set in file explorer in pickle format
    #     df.to_pickle("./temp_store/df_" + data_sample + "_" + str(nr_run) + ".pkl")

    # # Create and save multiple matrices
    # for sim_measure in sim_measures:
    #     all_matrices = []  # create list to store all matrices per measure
    #     for nr_run in range(total_nr_runs):
    #         df = pd.read_pickle("./temp_store/df_" + data_sample + "_" + str(nr_run) + ".pkl")
    #         matrix = create_sim_matrices(hpo_network, df, data_sample, sim_measure, nr_run)
    #         all_matrices.append(matrix)

    #     # Create and save average matrix per measure
    #     summed_matrix = np.zeros([len(df), len(df)])
    #     for m in all_matrices:
    #         summed_matrix = summed_matrix + m
    #     avg_matrix = summed_matrix / len(all_matrices)
    #     np.savetxt('./temp_store/sim_matrix_' + sim_measure + '_' + data_sample + '_avg', avg_matrix)

    # # OUTCOMMENT TILL HERE

    # Apply clustering and calculate ARI
    results = []
    for sim_measure in sim_measures:
        all_ARI = []
        for nr_run in range(total_nr_runs):
            df = pd.read_pickle("./temp_store/df_" + data_sample + "_" + str(nr_run) + ".pkl")
            ARI = cluster_sim_matrices(df, data_sample, sim_measure, nr_run, cluster_method)
            all_ARI.append(ARI)
        avg_ARI = sum(all_ARI)/len(all_ARI)
        stats_ARI = stats.describe(all_ARI)
        print("avg ARI: ", avg_ARI)
        # print("stats: ", stats_ARI)
        results.append([data_sample, sim_measure, avg_ARI, all_ARI, stats_ARI])

    for combi in measures_combi:
        # df = pd.read_pickle("./temp_store/df_" + data_sample + "_" + str(0) + ".pkl")
        # # this df does not correspond to created avg matrix, but the true
        # # labels are the same for all these data sets, which is what is used
        # ARI = cluster_sim_matrices(df, data_sample, combi, 'avg', cluster_method)
        # # run_nr is avg as the avg matrix needs to be used
        # results.append([data_sample, combi, ARI])

        all_ARI_combis = []
        for nr_run in range(total_nr_runs):
            df = pd.read_pickle("./temp_store/df_" + data_sample + "_" + str(nr_run) + ".pkl")
            ARI = cluster_sim_matrices(df, data_sample, combi, str(nr_run), cluster_method)
            all_ARI_combis.append(ARI)
        avg_ARI = sum(all_ARI_combis)/len(all_ARI_combis)
        stats_ARI = stats.describe(all_ARI_combis)
        print("avg ARI: ", avg_ARI)
        results.append([data_sample, combi, avg_ARI, all_ARI_combis, stats_ARI])

    return results


def create_sim_matrices(hpo_network, df, data_sample, sim_measure, nr_run):
    """
    Create similarity matrix according to a similarity measure and a data set.

    This function gets a certain data set and a string that indicates which
    similarity measure should be applied. The other parameters indicate the
    concerned data sample and the run number so that the matrix can be saved
    correctly.

    Parameters
    ----------
    hpo_network : classes.multidigraph.MultiDiGraph
        HPO network.
    df : DataFrame
        Dataframe that can be total data frame but also undersampled dataframe.
    data_sample : str
        String thatiindicates the concerned data set.
    sim_measure : str
        String that indicates the concerned similarity measure.
    nr_run : int
        Integer that indicates the concerned run number.

    Returns
    -------
    TYPE : Array of float64 or None
        Similarity matrix as a result of a given similarity measure and a
        certain data set. If given similarity measure is unknown return None.

    """
    # Create semantic similarity matrices
    if (sim_measure == 'Resnik'):
        print(sim_measure, " / Resnik matrix is created")
        scorer_Resnik = Scorer(hpo_network, scoring_method='Resnik')
        sim_matrix_Resnik = create_sem_sim_matrix(df, scorer_Resnik)
        np.savetxt('./temp_store/sim_matrix_' + sim_measure + '_' + data_sample + '_' + str(nr_run), sim_matrix_Resnik)
        return sim_matrix_Resnik

    if (sim_measure == 'Lin'):
        print(sim_measure, " / Lin matrix is created")
        scorer_Lin = Scorer(hpo_network, scoring_method='Lin')
        sim_matrix_Lin = create_sem_sim_matrix(df, scorer_Lin)
        np.savetxt('./temp_store/sim_matrix_' + sim_measure + '_' + data_sample + '_' + str(nr_run), sim_matrix_Lin)
        return sim_matrix_Lin

    # Create graph similarity matrices
    if (sim_measure == 'GED'):
        print(sim_measure, " / GED matrix is created")
        sim_matrix_GED = create_graph_sim_matrix_ged_gm(df)
        np.savetxt('./temp_store/sim_matrix_' + sim_measure + '_' + data_sample + '_' + str(nr_run), sim_matrix_GED)
        return sim_matrix_GED

    if (sim_measure == 'MCS'):
        print(sim_measure, " / MCS matrix is created")
        sim_matrix_MCS = create_graph_sim_matrix_mcs(df)
        np.savetxt('./temp_store/sim_matrix_' + sim_measure + '_' + data_sample + '_' + str(nr_run), sim_matrix_MCS)
        return sim_matrix_MCS

    else:
        print("No similarity measure was found. None instead of matrix is returned")
        return None


def cluster_sim_matrices(df, data_sample, sim_measure, nr_run, cluster_method):
    """
    Apply clustering algorithm to similarity matrices.

    Parameters
    ----------
    df : DataFrame
        Dataframe that can be total data frame but also undersampled dataframe.
    data_sample : str
        String thatiindicates the concerned data set.
    sim_measure : str
        String that indicates the concerned similarity measure.
    nr_run : int
        Integer that indicates the concerned run number.
    cluster_method : str
        String that indicates the which clustering method should be used.

    Returns
    -------
    ARI : float
        Adjusted Rand Index that indicates how good the clustering result is.
        Value between 0 and 1.

    """
    # Load all the similarity matrices

    # Cluster combined measures
    if (sim_measure == 'Resnik-MCS'):
        sim_matrix_Resnik = np.loadtxt(
            './temp_store/sim_matrix_Resnik_' + data_sample + '_' + str(nr_run))
        sim_matrix_MCS = np.loadtxt(
            './temp_store/sim_matrix_MCS_' + data_sample + '_' + str(nr_run))
        sim_matrix = (sim_matrix_Resnik + sim_matrix_MCS)/2

    elif (sim_measure == 'Lin-MCS'):
        sim_matrix_Lin = np.loadtxt(
            './temp_store/sim_matrix_Lin_' + data_sample + '_' + str(nr_run))
        sim_matrix_MCS = np.loadtxt(
            './temp_store/sim_matrix_MCS_' + data_sample + '_' + str(nr_run))
        sim_matrix = (sim_matrix_Lin + sim_matrix_MCS)/2

    elif (sim_measure == 'Resnik-Lin'):
        sim_matrix_Resnik = np.loadtxt(
            './temp_store/sim_matrix_Resnik_' + data_sample + '_' + str(nr_run))
        sim_matrix_Lin = np.loadtxt(
            './temp_store/sim_matrix_Lin_' + data_sample + '_' + str(nr_run))
        sim_matrix = (sim_matrix_Resnik + sim_matrix_Lin)/2

    elif (sim_measure == 'GED-MCS'):
        sim_matrix_GED = np.loadtxt(
            './temp_store/sim_matrix_GED_' + data_sample + '_' + str(nr_run))
        sim_matrix_MCS = np.loadtxt(
            './temp_store/sim_matrix_MCS_' + data_sample + '_' + str(nr_run))
        sim_matrix = (sim_matrix_GED + sim_matrix_MCS)/2

    # Cluster individual measures
    else:
        sim_matrix = np.loadtxt(
            './temp_store/sim_matrix_' + sim_measure + '_' + data_sample + '_' + str(nr_run))

    labels_true = df['label'].values
    # Transform into set to remove duplicates and take length:
    nr_labels = len(np.unique(labels_true))

    if (cluster_method == 'spectral'):
        ARI = cluster_kmeans(sim_measure, sim_matrix, labels_true, nr_labels)
    else:
        ARI = cluster_spectral(sim_measure, sim_matrix, labels_true, nr_labels)
    # show_similarity_matrix(measure_name, sim_matrix)
    return ARI


def calculate_variability(results, sim_measures, measures_combi):
    """
    Calculate variability of results of multiple applied runs.

    This function calculates the variability of the multiple ARIs obtained
    from the individual measures that have been applied multiple times to the
    data sets.

    Parameters
    ----------
    results : list
        List of lists including a.o. list of all ARIs, average ARI and stats.
    sim_measures : tuple
        Tuple with strings indicating the individual similarity measures.
    measures_combi : tuple
        Tuple with strings indicating the combined similarity measures.

    Returns
    -------
    variability_avg : float64
        The average variability calculated from minimum and maximum ARIs of
        every run.

    """
    var_sum = 0
    variabilities = []
    for r in results:
        # sum up the variability of every result by
        # subtracting the min from the max and adding it to the sum
        variabilities.append((r[4][1][0], r[4][1][1]))
        # variability_sum = var_sum + (r[4][1][1] - r[4][1][0])
    # divide total sum by length to get the avg
    # var_avg = var_sum/len(sim_measures)
    return variabilities


def calculate_pvals(results, sim_measures, measures_combi):
    """
    Calculate p-values of the differences between compared similarity measures.

    This function calculates p-values of the multiple ARIs obtained from the
    measures that have been applied multiple times to the data sets.
    Therefore, all the values of one measure are compared with all the values
    of the other measure, to which the stats.ttest is applied.

    Parameters
    ----------
    results : list
        List of lists including a.o. list of all ARIs, average ARI and stats.
    sim_measures : tuple
        Tuple with strings indicating the individual similarity measures.
    measures_combi : tuple
        Tuple with strings indicating the combined similarity measures.

    Returns
    -------
    pvals : list
        List of tuples that consist of the two measures compared and the
        resulting p-value that indicates if the difference is significant.

    """
    pvals = []
    for i, r1 in enumerate(results):
        if not(any(r1[1] == measure[0] for measure in pvals)):
            for j, r2 in enumerate(results):
                if ((not(r2 == r1)) and (not(any(r2[1] == measure[0]
                                                 for measure in pvals)))):
                    statist, pval = stats.ttest_ind(r1[3], r2[3])
                    pvals.append((r1[1], r2[1], pval))
    return pvals


# def main():

df_updated_total, hpo_network = prepare_data()

# # Save to and get dataframe from pickle
# df_updated_total.to_pickle("data_total_updated.pkl")
# df_updated_total = pd.read_pickle("data_total_updated.pkl")

data_samples = (
                ('total'),
                ('undersampled'),
                ('undersampled_withoutspop'),
                )

sim_measures = (
                ('Resnik'),
                ('Lin'),
                ('GED'),
                ('MCS'),
                )

measures_combi = (
                ('Resnik-MCS'),  # combine best semantic and graph measure
                ('Lin-MCS'),
                ('Resnik-Lin'),  # combine both semantic measures
                ('GED-MCS'),  # combine both graph measures
                )

cluster_method = 'kmeans'  # indicate cluster method here
total_nr_runs = 20  # indicate number of runs for undersampled data here

# Create empty lists to store values
results_total_dataset = []
results_mtp_runs_undersampled = []
results_mtp_runs_undersampled_withoutspop = []
variability_undersampled = []
variability_undersampled_withoutspop = []
pvals_undersampled = []
pvals_undersampled_withoutspop = []

# Loop through the diffent data sets and apply the whole method including
# all the different (combinations of) similarity measures to every data
# set. For the undersampled data sets multiple runs need to be executed.
for data_sample in data_samples:
    if (data_sample == 'total'):
        results_total_dataset = pipeline_single_run(
            df_updated_total,
            hpo_network,
            data_sample,
            sim_measures,
            measures_combi,
            cluster_method)

    elif(data_sample == 'undersampled'):
        results_mtp_runs_undersampled = pipeline_mtp_runs(
            df_updated_total,
            hpo_network,
            data_sample,
            sim_measures,
            measures_combi,
            cluster_method,
            total_nr_runs)

        variability_undersampled = calculate_variability(
            results_mtp_runs_undersampled, sim_measures, measures_combi)
        pvals_undersampled = calculate_pvals(
            results_mtp_runs_undersampled, sim_measures, measures_combi)

    elif(data_sample == 'undersampled_withoutspop'):
        results_mtp_runs_undersampled_withoutspop = pipeline_mtp_runs(
            df_updated_total,
            hpo_network,
            data_sample,
            sim_measures,
            measures_combi,
            cluster_method,
            total_nr_runs)

        variability_undersampled_withoutspop = calculate_variability(
            results_mtp_runs_undersampled_withoutspop, sim_measures, measures_combi)
        pvals_undersampled_withoutspop = calculate_pvals(
            results_mtp_runs_undersampled_withoutspop, sim_measures, measures_combi)

    else:
        print("The", data_sample, "data sample is unknown.")


# if __name__ == "__main__":
#     main()
