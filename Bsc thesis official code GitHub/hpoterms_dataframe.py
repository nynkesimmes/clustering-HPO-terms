# -*- coding: utf-8 -*-
"""
Created on Wed Nov 17 11:12:40 2021

@author: nynkesimmes
"""

def update_df(df, obo_file_path):
    """
    Update original dataframe by calling helper functions.

    The function checks for duplicates in each list of HPO terms.
    A dictionary is created of all terms that have been replaced with other
    HPO id's. Then is looped through all the terms and checked by use of the
    dictionary if their id needs to be replaced. In the end all
    the non_phenotypes are removed.

    Parameters
    ----------
    df : DataFrame
        Dataframe of 553 patients with HPO terms and corresponding HPO graphs.
    obo_file_path : str
        Path to hpo.obo file.

    Returns
    -------
    df_updated : DataFrame
        Completely updated dataframe.

    """

    df_without_dupl = remove_dupl_terms(df, obo_file_path)  # remove duplicates

    # Create dictionary and replace old hpo terms with their alternative ids
    df_updated = df_without_dupl.copy()
    id_dict = create_dict_replaced_ids(obo_file_path)
    print("Replacing old hpo terms with their alternative id")
    for i, terms in enumerate(df_without_dupl['hpo_all']):
        for j, term in enumerate(terms):
            if term in id_dict.keys():
                df_updated['hpo_all'][i][j] = id_dict.get(term)
            if (term == 'HP:0006830'):  # hardcoded, is only alternative id
                df_updated['hpo_all'][i][j] = 'HP:0001319'

    # Remove all hpo terms that are non_phenotypes; these are also removed
    # when the network is generated
    df_updated2 = remove_non_phenotypes(df_updated, obo_file_path)
    return df_updated2


def remove_dupl_terms(df, obo_file_path):
    """
    Remove duplicate HPO terms in each list of HPO terms of patients.

    Parameters
    ----------
    df : DataFrame
        Dataframe of 553 patients with HPO terms and corresponding HPO graphs.
    obo_file_path : str
        Path to hpo.obo file.

    Returns
    -------
    df_without_dupl : DataFrame
        Dataframe without duplicates in each list of HPO terms.

    """
    print("Removing duplicate hpo terms")
    df_without_dupl = df.copy()
    for i, terms in enumerate(df['hpo_all']):
        terms_set = set(terms)
        terms_without_dupl = list(terms_set)  # duplicates a.o.: 77, 88, 90
        df_without_dupl['hpo_all'][i] = terms_without_dupl
    return df_without_dupl


def create_dict_replaced_ids(obo_file_path):
    """
    Create dictionary with all to id's that are replaced by a new id.

    In this dictionary all keys represent the "old" id values and the values
    represent the id's they should be replaced with.

    Parameters
    ----------
    obo_file_path : str
        Path to hpo.obo file.

    Returns
    -------
    id_dict : dict
        Dictionary of all to be replaced id's and the new id's.

    """
    print("Create dictionary with ids of hpo terms that need to be replaced")
    id_dict = {}
    for line in open(obo_file_path):
        if line.startswith("id: "):
            current_id = line[len("id: "):-1]  # get the part of the line after "id:" and until "\n" (therefore -1)
        if line.startswith("replaced_by"):
            new_id = line[len("replaced_by: "):-1]
            id_dict[current_id] = new_id
    return id_dict


def remove_non_phenotypes(df_updated, obo_file_path):
    """
    Remove non_phenotypes and its child nodes from dataframe.

    Parameters
    ----------
    df : DataFrame
        Dataframe of 553 patients with HPO terms and corresponding HPO graphs.
    obo_file_path : str
        Path to hpo.obo file.

    Returns
    -------
    df_without_dupl : DataFrame
        Dataframe without non_phenotypes and its child nodes.
    """
    print("Removing all hpo terms that are non-phenotypes")
    df_updated2 = df_updated.copy()
    all_non_phenotypes = get_all_non_phenotypes(obo_file_path)

    for i, terms in enumerate(df_updated['hpo_all']):
        for j, term in enumerate(terms):
            if term in all_non_phenotypes:
                # print("Non-phenotype: ", term)
                del df_updated2['hpo_all'][i][j]
    return df_updated2


def get_all_non_phenotypes(obo_file_path):
    """
    Create list of all non_phenotypes and its child nodes.

    Parameters
    ----------
    obo_file_path : str
        Path to hpo.obo file.

    Returns
    -------
    non_phenotypes : list
        List of all non_phenotypes and its child nodes.

    """
    # roots for non-phenotype nodes
    non_phenotypes = ['HP:0040006', 'HP:0000005', 'HP:0012823', 'HP:0040279',
                      'HP:0031797']
    nr_ids_added = 1
    while(nr_ids_added > 0):  # while new ids are added
        nr_ids_added = 0
        for line in open(obo_file_path):
            if line.startswith("id: "):
                current_id = line[len("id: "):-1]
            if line.startswith("is_a: "):
                is_a_term = line[len("is_a: "):len("is_a: ")+10]  # 10 is length of HPO term
                if ((is_a_term in non_phenotypes) & (current_id not in non_phenotypes)):
                    non_phenotypes.append(current_id)
                    nr_ids_added += 1
    return non_phenotypes
