#!/usr/bin/env python

"""
tcrd.py: Query Target Central Resource Database (TCRD).
"""

__author__ = "Maria J. Falaguera"
__date__ = "21 Oct 2024"

import pymysql
import pandas as pd

# from rdkit import Chem


dbDir = None  #'/mnt/ct-server/DB/sources/affinityDBs/ChemblDB/2021-28'

cutoffs = {
    "Kinase": 7.7,
    "GPCR": 6.5,
    "oGPCR": 6.5,
    "NR": 6.1,
    "Transporter": 6.1,
    "Enzyme": 5.2,
    "Cytochrome P450": 5.2,
    "IC": 4.6,
    "Other": 6.3,
    "Epigenetic": 6.3,
    "TF; Epigenetic": 6.3,
    "TF": 6.3,
    "NAN": 100,
}

names = [
    ["Enzyme", "Enzymes"],
    ["Epigenetic", "Epig. proteins"],
    ["GPCR", "GPCRs"],
    ["IC", "Ion channels"],
    ["Kinase", "Kinases"],
    ["NR", "Nuc. rec."],
    ["TF", "Transc. factors"],
    ["Transporter", "Transporters"],
    ["Other", "Others"],
    ["inter-family", "Inter-family"],
]

tdls = ["Tclin", "Tchem"]

order = [[name, idx] for idx, [_, name] in enumerate(names)]
tdl_order = [[name, idx] for idx, name in enumerate(tdls)]

palette = {
    "TTclin": "#092957",
    "TTchem": "#1A9849",
    "TTmix": "#FD9407",
    "DID(MOA)": "#f6b81c",
    "DID(MIX)": "#ca387d",
    "DID(OFF)": "#207bbd",
    "DIL(LIG)": "#c9c9c9",
    "Tclin on-target": "#092957",
    "Tclin off-target": "#508feb",
    "Tchem off-target": "#1A9849",
    "[0, 1)": "#207bbd",
    "[1, 2)": "#ca387d",
    "≥2": "#f6b81c",
}

selectivities = {0: "[0, 1)", 1: "[1, 2)", 2: "≥2"}


def openConnection():
    # according to slack chat with Eloy, this tcrd version must be:
    #  TCRD v6.12.4
    # Release Date: 20211029
    conn = pymysql.connect(
        database="tcrd",
        user="tcrd",
        host="mysql-chembl-mychembl-dev.ebi.ac.uk",
        password="tcrd",
        port=4233,
    )
    return conn


def getProteins():
    """
    Get proteins.
    """

    SELECT = ["protein.id", "protein.uniprot"]
    FROM = ["protein"]
    WHERE = ["1 = 1"]

    sql = "SELECT DISTINCT {} FROM {} WHERE {};".format(
        ",".join(SELECT), ",".join(FROM), " AND ".join(WHERE)
    )

    results = {}
    conn = openConnection()
    cur = conn.cursor()
    cur.execute(sql)
    for idx, uniprot in cur:
        results[idx] = uniprot

    cur.close()
    conn.close()

    return results


def countProteins():
    return len(getProteins())


def getTargets():
    """
    Get targets.
    """

    SELECT = ["target.id"]
    FROM = ["target"]
    WHERE = ["1 = 1"]

    sql = "SELECT DISTINCT {} FROM {} WHERE {};".format(
        ",".join(SELECT), ",".join(FROM), " AND ".join(WHERE)
    )

    results = set()
    conn = openConnection()
    cur = conn.cursor()
    cur.execute(sql)
    for idx in cur:
        results.update([idx[0]])

    cur.close()
    conn.close()

    return results


def getInfoForProtein(uniprot=None, output=None):
    """
    Get info for protein.

    Args:
        uniprot (str/list/set/dict):    protein uniprot id
        idx (str/list/set/dict):        protein ID in TCRD

    Returns:
        dict:   {uniprot: {'tdl': tdl, 'name': name, 'family': family, ...}}
    """

    SELECT = [
        "p.uniprot as targetUniprot",
        "t.tdl as targetTDL",
        "t.fam as targetFamily",
        "p.sym as targetSymbol",
        "p.description as targetName",
        # "p.seq as targetSequence",
        # "p.stringid as gene_id",
        "p.family as targetSubfamily",
        # "t.ttype",
        # "p.dtoid",
        # "p.dtoclass",
        # "p.id as protein_id",
        "t.id as targetId",
    ]

    FROM = ["protein p"]
    INNER_JOIN = ["t2tc ON t2tc.protein_id = p.id", "target t ON t2tc.target_id = t.id"]
    WHERE = ["1 = 1"]

    if uniprot is not None:
        if isinstance(uniprot, str):
            uniprot = [uniprot]
        WHERE.append("p.uniprot IN ('{}')".format("','".join(list(uniprot))))

    sql = "SELECT DISTINCT {} FROM {} INNER JOIN {} WHERE {};".format(
        ",".join(SELECT),
        ",".join(FROM),
        " INNER JOIN ".join(INNER_JOIN),
        " AND ".join(WHERE),
    )

    results = {}
    conn = openConnection()
    cur = conn.cursor(pymysql.cursors.DictCursor)
    cur.execute(sql)
    results = [{r: row[r] for r in row} for row in cur]

    cur.close()
    conn.close()

    if output == "pandas":
        results = pd.DataFrame(results)

    return results


def getTDLforProtein(uniprot=None, by="gene_name"):
    """
    Get TDL for protein.

    Args:
        uniprot (str/list/set/dict):    protein uniprot id

    Returns:
        dict:   {uniprot: TDL}
    """

    results = {r[by]: r["tdl"] for r in getInfoForProtein(uniprot=uniprot)}
    return results


def getNameForProtein(uniprot=None):
    """
    Get name for protein.

    Args:
        uniprot (str/list/set/dict):    protein uniprot id

    Returns:
        dict:   {uniprot: name}
    """

    results = {
        t: info["name"] for t, info in getInfoForProtein(uniprot=uniprot).items()
    }
    return results


def getDescriptionForProtein(uniprot=None):
    """
    Get description for protein.

    Args:
        uniprot (str/list/set/dict):    protein uniprot id

    Returns:
        dict:   {uniprot: description}
    """

    results = {
        t: info["description"] for t, info in getInfoForProtein(uniprot=uniprot).items()
    }
    return results


def getGeneForProtein(uniprot=None):
    """
    Get gene for protein.

    Args:
        uniprot (str/list/set/dict):    protein uniprot id

    Returns:
        dict:   {uniprot: gene}
    """

    results = {
        t: info["gene_name"] for t, info in getInfoForProtein(uniprot=uniprot).items()
    }
    return results


def getProteinFamilyForProtein(uniprot=None):
    """
    Get protein family for protein.

    Args:
        uniprot (str/list/set/dict):    protein uniprot id

    Returns:
        dict:   {uniprot: protein_family}
    """

    results = {
        t: info["protein_family"]
        for t, info in getInfoForProtein(uniprot=uniprot).items()
    }
    return results


def getChemblRefYear(uniprot=None, include_empty=False, by="gene_symbol"):
    """
    Get ChEMBL First Reference Year for protein.

    Args:
        uniprot (str/list/set/dict):    protein uniprot id
        include_empty (bool):           include or don't include in the return proteints w/o ChEMBL First Reference Year

    Returns:
        dict:   {uniprot: chembl_ref_year, ...}}
    """

    SELECT = [
        "p.uniprot",
        "p.sym as gene_symbol",
        "tdl_info.integer_value as chembl_ref_year",
    ]
    FROM = ["protein p"]
    LEFT_JOIN = [
        "t2tc ON t2tc.protein_id = p.id",
        "tdl_info ON tdl_info.target_id = t2tc.target_id",
    ]
    WHERE = ['tdl_info.itype = "ChEMBL First Reference Year"']

    if uniprot is not None:
        if isinstance(uniprot, str):
            uniprot = [uniprot]
        WHERE.append("p.uniprot IN ('{}')".format("','".join(list(uniprot))))

    sql = "SELECT DISTINCT {} FROM {} LEFT JOIN {} WHERE {};".format(
        ",".join(SELECT),
        ",".join(FROM),
        " INNER JOIN ".join(LEFT_JOIN),
        " AND ".join(WHERE),
    )

    results = {}
    conn = openConnection()
    cur = conn.cursor(pymysql.cursors.DictCursor)
    cur.execute(sql)
    for row in cur:
        results[row[by]] = row["chembl_ref_year"]

    cur.close()
    conn.close()

    if include_empty:
        for p in set(uniprot).difference(set(results)):
            results[p] = None

    return results


def getChemblSelectiveCompound(uniprot=None):
    """
    Get ChEMBL Selective Compound for protein.

    Args:
        uniprot (str/list/set/dict):    protein uniprot id

    Returns:
        dict:   {uniprot: chembl_sel_cmpd, ...}}
    """

    SELECT = ["p.uniprot", "tdl_info.string_value as chembl_sel_cmpd"]
    FROM = ["protein p"]
    INNER_JOIN = [
        "t2tc ON t2tc.protein_id = p.id",
        "tdl_info ON tdl_info.target_id = t2tc.target_id",
    ]
    WHERE = ['tdl_info.itype = "ChEMBL Selective Compound"']

    if uniprot is not None:
        if isinstance(uniprot, str):
            uniprot = [uniprot]
        WHERE.append("p.uniprot IN ('{}')".format("','".join(list(uniprot))))

    sql = "SELECT DISTINCT {} FROM {} INNER JOIN {} WHERE {};".format(
        ",".join(SELECT),
        ",".join(FROM),
        " INNER JOIN ".join(INNER_JOIN),
        " AND ".join(WHERE),
    )

    results = {}
    conn = openConnection()
    cur = conn.cursor(pymysql.cursors.DictCursor)
    cur.execute(sql)
    for row in cur:
        results[row["uniprot"]] = row["chembl_sel_cmpd"]

    cur.close()
    conn.close()

    return results


def getTargetFamilyForProtein(uniprot=None, agglutinate=False):
    """
    Get family for protein.

    Args:
        uniprot (str/list/set/dict):    protein uniprot id
        agglutinate (bbol):             convert oGPCR into GPCR, TF; Epigenetic into Epigenetic

    Returns:
        dict:   {uniprot: family}
    """

    results = {
        t: info["target_family"]
        for t, info in getInfoForProtein(uniprot=uniprot).items()
    }

    if agglutinate:
        results2 = {}
        for t in results:
            family = results[t]
            if family == "oGPCR":
                family = "GPCR"
            elif family == "TF; Epigenetic":
                family = "Epigenetic"
            results2[t] = family
        results = results2

    return results


def getProteinsForFamily(family=None):
    """
    Get proteins for protein family.

    Args:
        family (str/list/set/dict):    protein family

    Returns:
        dict:   {family: set of proteins, ...}
    """

    SELECT = ["family", "uniprot"]
    FROM = ["protein"]

    if family is None:
        results = {}
        conn = openConnection()
        sql = "SELECT DISTINCT {} FROM {};".format(",".join(SELECT), ",".join(FROM))
        cur = conn.cursor(pymysql.cursors.DictCursor)
        cur.execute(sql)
        for row in cur:
            try:
                results[row["family"]].update([row["uniprot"]])
            except KeyError:
                results[row["family"]] = set([row["uniprot"]])
        cur.close()
        conn.close()

    elif isinstance(family, str):
        results = {}
        conn = openConnection()
        sql = 'SELECT DISTINCT {} FROM {} WHERE family = "{}";'.format(
            ",".join(SELECT), ",".join(FROM), family
        )
        cur = conn.cursor(pymysql.cursors.DictCursor)
        cur.execute(sql)
        for row in cur:
            try:
                results[row["family"]].update([row["uniprot"]])
            except KeyError:
                results[row["family"]] = set([row["uniprot"]])
        cur.close()
        conn.close()

    else:
        results = {}
        conn = openConnection()
        i = 0
        n = 5

        while i < len(family):
            sub_family = family[i : i + n]
            subfamily = "("
            for f in family[i : i + n]:
                if '"' in f:
                    subfamily += "'{}',".format(f)
                elif "'" in f:
                    subfamily += '"{}",'.format(f)
                else:
                    subfamily += '"{}",'.format(f)
            subfamily = subfamily.strip(",") + ")"

            WHERE = ["family IN {}".format(subfamily)]

            sql = "SELECT DISTINCT {} FROM {} WHERE {};".format(
                ",".join(SELECT), ",".join(FROM), " AND ".join(WHERE)
            )

            cur = conn.cursor(pymysql.cursors.DictCursor)
            cur.execute(sql)
            for row in cur:
                try:
                    results[row["family"]].update([row["uniprot"]])
                except KeyError:
                    results[row["family"]] = set([row["uniprot"]])
            cur.close()

            i += n

        conn.close()

    return results


def getSiblingsForProtein(uniprot=None):
    """
    Get siblings (within the same protein family) for protein.

    Args:
        uniprot (str/list/set/dict):    protein uniprot id

    Returns:
        dict:   {uniprot: set of sibling uniprots}
    """

    protein_family = getProteinFamilyForProtein(uniprot=uniprot)
    families = set(protein_family.values()).difference(set([None]))
    family_proteins = getProteinsForFamily(family=families)

    results = {}  # protein -> siblings
    for protein in protein_family:
        family = protein_family[protein]
        if family is None:
            siblings = None
        else:
            siblings = family_proteins[family].difference(set([protein]))
        results[protein] = siblings

    return results


def getCmpdActivitiesForProtein(
    uniprot=None,
    smiles=False,
    pact_cutoff=False,
):
    """
    Get compound activities for protein.

    Args:
        uniprot (str/list/set/tuple/dict):    protein to query
        smiles (bool):                        only activities implying cmpd with SMILES
        pact_cutoff (boold/float/int):        only activities with a value above cutoff
        pubmed_year (bool):                   annotate pubmed_ids with their publication year

    Returns:
        dict: { uniprot: {'cmpd_id': [{'act_value': act_value, 'pubmed_ids': pubmed_ids set}, ...], ... }
    """

    WHERE = ["1 = 1"]

    if smiles:
        WHERE.append("cmpd_activity.smiles != ''")

    if pact_cutoff:
        WHERE.append("cmpd_activity.act_value >= {}".format(pact_cutoff))

    if uniprot is not None:
        if isinstance(uniprot, str):  # single target
            WHERE.append("p.uniprot = '{}'".format(uniprot))
        elif (
            isinstance(uniprot, set)
            or isinstance(uniprot, list)
            or isinstance(uniprot, tuple)
            or isinstance(uniprot, dict)
        ):  # multiple targets
            WHERE.append("p.uniprot IN ('{}')".format("','".join(list(uniprot))))

    # Get cmpd activities
    SELECT = [
        "cmpd_activity.target_id",
        "cmpd_activity.cmpd_pubchem_cid",
        "cmpd_activity.act_value",
        "cmpd_activity.act_type",
        # "cmpd_activity.pubmed_ids",
    ]
    FROM = ["cmpd_activity"]

    sql = "SELECT {} FROM {} WHERE {};".format(
        ",".join(SELECT), ",".join(FROM), " AND ".join(WHERE)
    )

    conn = openConnection()
    cur = conn.cursor(pymysql.cursors.DictCursor)
    cur.execute(sql)

    results = [{r: row[r] for r in row} for row in cur]

    cur.close()
    conn.close()

    return results


def getDrugActivitiesForProtein(
    uniprot=None, smiles=False, pact_cutoff=False, has_moa=False
):
    """
    Get compound activities for protein.

    Args:
        uniprot (str/list/set/tuple/dict):    protein to query
        smiles (bool):                        only activities implying cmpd with SMILES
        pact_cutoff (boold/float/int):        only activities with a value above cutoff

    Returns:
        dict: { uniprot: {'cmpd_id': [{'act_value': act_value, 'pubmed_ids': pubmed_ids set}, ...], ... }
    """

    WHERE = ["1 = 1"]

    if smiles:
        WHERE.append("drug_activity.smiles IS NOT NULL")

    if pact_cutoff:
        WHERE.append("drug_activity.act_value >= {}".format(pact_cutoff))

    if has_moa == True:
        WHERE.append("drug_activity.has_moa = 1")

    if uniprot is not None:
        if isinstance(uniprot, str):  # single target
            WHERE.append("protein.uniprot = '{}'".format(uniprot))
        elif (
            isinstance(uniprot, set)
            or isinstance(uniprot, list)
            or isinstance(uniprot, tuple)
            or isinstance(uniprot, dict)
        ):  # multiple targets
            WHERE.append("protein.uniprot IN ('{}')".format("','".join(list(uniprot))))

    # Get drug activities
    SELECT = [
        "drug_activity.target_id",
        "drug_activity.cmpd_pubchem_cid",
        "drug_activity.act_value",
        "drug_activity.act_type",
        "drug_activity.has_moa",
        # "disease.name AS indication",
    ]
    FROM = ["drug_activity"]
    LEFT_JOIN = [
        "t2tc ON t2tc.target_id = drug_activity.target_id",
        "disease ON disease.protein_id = t2tc.protein_id AND disease.drug_name = drug_activity.drug",
        "protein ON protein.id = t2tc.protein_id",
    ]

    sql = "SELECT {} FROM {} LEFT JOIN {} WHERE {};".format(
        ",".join(SELECT),
        ",".join(FROM),
        " LEFT JOIN ".join(LEFT_JOIN),
        " AND ".join(WHERE),
    )
    print(sql)

    conn = openConnection()
    cur = conn.cursor(pymysql.cursors.DictCursor)
    cur.execute(sql)

    results = [{r: row[r] for r in row} for row in cur]

    cur.close()

    return results


def getCutoffForFamily(family=None):
    """
    Dictionary of family-specific pAct cutoffs.

    Args:
        family (str/set/tuple/dict, default=None):  Family to get cutoff

    Returns:
        dict:   Dictionary with families as keys and cutoffs as values
    """

    cutoffs = {
        "Kinase": 7.5,  # <= 30nM
        "GPCR": 7.0,  # <= 100nM
        "IC": 5.0,  # <= 10uM
        "Other": 6.0,  # <= 1uM (non-IDG family)
    }

    cutoffs["KC"] = cutoffs["Kinase"]
    cutoffs["GR"] = cutoffs["GPCR"]
    cutoffs["Ion channel"] = cutoffs["IC"]

    if family is not None:
        if family not in cutoffs:
            cutoffs[family] = cutoffs["Other"]

    return cutoffs


def getColorForTDL():
    """
    Dictionary of TDL-specific colors.
    """

    colors = {
        "Tclin": "#01437f",
        "Tchem": "#23b061",
        "Tbio": "#f26a6c",
        "Tdark": "#242021",
    }

    return colors


def getColorForFamily(alpha=1):
    """
    Dictionary of family-specific colors.

    Args:
        alpha (float):    transparency

    Returns:
        dict:   Dictionary with families as keys and colors as values
    """

    colors = {
        "GPCR": "#8eb566",
        "Nuclear receptor": "#b00b29",
        "Ion channel": "#248ec8",
        "Enzyme": "#8e84b6",
        "Transporter": "#f0ef63",
        "Kinase": "#d65192",
        "Transcription factor": "#ea7e1c",
        "Epigenetic": "#86c8d8",
        "Other": "#534b49",
    }

    colors = {
        "GPCR": (200, 235, 164, alpha),
        "Nuclear receptor": (176, 11, 41, alpha),
        "Ion channel": (36, 142, 200, alpha),
        "Enzyme": (142, 132, 182, alpha),
        "Transporter": (240, 239, 99, alpha),
        "Kinase": (214, 81, 146, alpha),
        "Transcription factor": (234, 126, 28, alpha),
        "Epigenetic": (134, 200, 216, alpha),
        "Other": (83, 75, 73, alpha),
    }

    light_colors = {
        "GPCR": "#d8f0c0",  #'#c0dea2',
        "Nuclear receptor": "#fcc7d1",  ##ffa6b6',
        "Ion channel": "#d9f1ff",  #'#b3e4ff',
        "Enzyme": "#eae3ff",  ##d6c9ff',
        "Transporter": "#ffffba",  #'#f0ef63',
        "Kinase": "#ffd4e9",  ##ffabd4',
        "Transcription factor": "#ffd6b0",  ##ea7e1c',
        "Epigenetic": "#dbf8ff",  ##c2f3ff',
        "Other": "#e6cfca",  ##534b49'
    }

    colors["GR"] = colors["GPCR"]
    colors["NR"] = colors["Nuclear receptor"]
    colors["IC"] = colors["Ion channel"]
    colors["EC"] = colors["Enzyme"]
    colors["TR"] = colors["Transporter"]
    colors["KC"] = colors["Kinase"]
    colors["TF"] = colors["Transcription factor"]
    colors["EP"] = colors["Epigenetic"]

    return colors


def getInfoForPubmed(pubmed_id=None):
    """
    Get info for pubmed ID.

    Args:
        pubmed_id (str/list/set/dict):    pubmed_if to query

    Returns:
        dict:   {pubmed_id: {'year': year, ...}, ...}
    """

    SELECT = ["*"]
    FROM = ["pubmed"]
    WHERE = ["1 = 1"]

    if pubmed_id is not None:
        if isinstance(pubmed_id, int):
            WHERE.append("id = {}".format(pubmed_id))
        elif (
            isinstance(pubmed_id, list)
            or isinstance(pubmed_id, tuple)
            or isinstance(pubmed_id, dict)
            or isinstance(pubmed_id, set)
        ):
            if len(pubmed_id):
                WHERE.append("id IN ({})".format(",".join([str(p) for p in pubmed_id])))
            else:  # empty list passed from getCmpdActivityForProtein
                return {}

    sql = "SELECT DISTINCT {} FROM {} WHERE {};".format(
        ",".join(SELECT), ",".join(FROM), " AND ".join(WHERE)
    )

    results = {}
    conn = openConnection()
    cur = conn.cursor(pymysql.cursors.DictCursor)
    cur.execute(sql)
    for row in cur:
        results[row["id"]] = {r: row[r] for r in row if r != "id"}

        # add publication year
        if (row["date"] is None) or (row["date"] == "None"):
            year = None
        else:
            year = int(row["date"].split("-")[0])
        results[row["id"]]["year"] = year

    cur.close()
    conn.close()

    return results


def getYearForPubmed(pubmed_id=None):
    """
    Get year for pubmed_id.

    Args:
        pubmed_id (str/list/set/dict):    pubmed_id to query

    Returns:
        dict:   {pubmed_id: year}
    """

    results = {
        p: info["year"] for p, info in getInfoForPubmed(pubmed_id=pubmed_id).items()
    }
    return results


def getInfoForCmpd(cmpd=None):
    """
    Get info for cmpd.

    Args:
        cmpd (str/list/set/dict):    cmpd_id_in_src to query

    Returns:
        dict:   {cmpd_id_in_src: {'tdl': tdl, 'name': name, 'family': family, ...}}
    """

    SELECT = ["cmpd_name_in_src", "smiles"]
    FROM = ["cmpd_activity"]
    WHERE = ["1 = 1"]

    if cmpd is not None:
        if isinstance(cmpd, str):
            cmpd = [cmpd]
        WHERE.append("cmpd_id_in_src IN ('{}')".format("','".join(list(cmpd))))

    sql = "SELECT DISTINCT {} FROM {} WHERE {};".format(
        ",".join(SELECT), ",".join(FROM), " AND ".join(WHERE)
    )

    results = {}
    conn = openConnection()
    cur = conn.cursor(pymysql.cursors.DictCursor)
    cur.execute(sql)
    results = [{r: row[r] for r in row} for row in cur]

    cur.close()
    conn.close()

    return results


def getSmilesForCmpd(cmpd=None):
    """
    Get smiles for cmpd.

    Args:
        cmpd (str/list/set/dict):    cmpd_id_in_src to query

    Returns:
        dict:   {cmpd_id_in_src to query: smiles}
    """

    results = [
        {
            "smiles": r["smiles"],
            "cmpd_id_in_src": r["cmpd_id_in_src"],
            "cmpd_pubchem_cid": r["cmpd_pubchem_cid"],
        }
        for r in getInfoForCmpd(cmpd=cmpd)
    ]
    return results


def getInchikeyForCmpd(cmpd=None):
    """
    Get smiles for cmpd.

    Args:
        cmpd (str/list/set/dict):    cmpd_id_in_src to query

    Returns:
        dict:   {cmpd_id_in_src to query: inchikey}
    """

    results = {}
    for c, info in getInfoForCmpd(cmpd=cmpd).items():
        mol = Chem.MolFromSmiles(info["smiles"])
        if mol is not None:
            inchi = Chem.inchi.MolToInchi(mol)
            inchikey = Chem.inchi.InchiToInchiKey(inchi)
        else:
            print(info["smiles"])
            inchikey = None
        results[c] = inchikey

    return results


def getInfoForDrug(drug=None):
    """
    Get info for drug.

    Args:
        drug (str/list/set/dict):    drug to query

    Returns:
        dict:   {drug: {'smiles': smiles}, ...}
    """

    SELECT = ["drug", "smiles"]
    FROM = ["drug_activity"]
    WHERE = ["1 = 1"]

    if drug is not None:
        if isinstance(drug, str):
            drug = [drug]
        WHERE.append("drug IN ('{}')".format("','".join(list(drug))))

    sql = "SELECT DISTINCT {} FROM {} WHERE {};".format(
        ",".join(SELECT), ",".join(FROM), " AND ".join(WHERE)
    )

    results = {}
    conn = openConnection()
    cur = conn.cursor(pymysql.cursors.DictCursor)
    cur.execute(sql)
    results = [{r: row[r] for r in row} for row in cur]

    cur.close()
    conn.close()

    return results


def getSmilesForDrug(drug=None):
    """
    Get smiles for drug.

    Args:
        drug (str/list/set/dict):    drug to query

    Returns:
        dict:   {drug: smiles, ...}
    """

    results = [
        {
            "smiles": r["smiles"],
            "drug": r["drug"],
            "cmpd_pubchem_cid": r["cmpd_pubchem_cid"],
        }
        for r in getInfoForDrug(drug=drug)
    ]
    return results


def getInchikeyForDrug(drug=None):
    """
    Get inchikey for drug.

    Args:
        drug (str/list/set/dict):    drug to query

    Returns:
        dict:   {drug: inchikey, ...}
    """

    results = {}
    for d, info in getInfoForDrug(drug=drug).items():
        mol = Chem.MolFromSmiles(info["smiles"])
        if mol is not None:
            inchi = Chem.inchi.MolToInchi(mol)
            inchikey = Chem.inchi.InchiToInchiKey(inchi)
        else:
            inchikey = None
        results[d] = inchikey

    return results


def getIndicationForTarget():
    SELECT = [
        "disease.name as indication",
        "t2tc.target_id",
    ]
    FROM = ["disease"]
    LEFT_JOIN = ["t2tc ON t2tc.protein_id = disease.protein_id"]

    sql = "SELECT DISTINCT {} FROM {} LEFT JOIN {};".format(
        ",".join(SELECT), ",".join(FROM), " LEFT JOIN ".join(LEFT_JOIN)
    )

    conn = openConnection()
    cur = conn.cursor(pymysql.cursors.DictCursor)
    cur.execute(sql)

    results = [{r: row[r] for r in row} for row in cur]

    cur.close()
    conn.close()

    return results


## ----- NOT USED IN ILLUMINATING THE CHEMCIAL SPACE (ICS) PROJECT -------- ##


def filterActivitiesForProtein(
    activities, by_act_cutoff=False, by_year_threshold=False
):
    tmp = {}

    # by_act_cutof is not False
    if by_act_cutoff:
        # if by_act_cutoff is True -> filter by IDG-family cutoff
        if isinstance(by_act_cutoff, bool):
            protein_family = getTargetFamilyForProtein(uniprot=set(activities))
            family_cutoff = getCutoffForFamily(family=None)

            for protein in activities:
                # get family cutoff
                family = protein_family[protein]
                try:
                    cutoff = family_cutoff[family]
                except KeyError:
                    cutoff = family_cutoff["Other"]  # non-IDG family

                # print(protein, family, cutoff)

                for cmpd in activities[protein]:
                    for data in activities[protein][cmpd]:
                        act_value = data["act_value"]
                        pubmed_ids = data["pubmed_ids"]

                        if act_value >= cutoff:
                            # add protein
                            try:
                                _ = tmp[protein]
                            except KeyError:  # new protein
                                tmp[protein] = {}

                            # add cmpd for protein
                            try:
                                _ = tmp[protein][cmpd]
                            except KeyError:  # new cmpd for protein
                                tmp[protein][cmpd] = []

                            # add pubmed_ids and act_value for cmpd for protein
                            tmp[protein][cmpd].append(
                                {"act_value": act_value, "pubmed_ids": pubmed_ids}
                            )

        # if by_act_cutoff is a value -> filter by this act cutoff
        else:
            cutoff = by_act_cutoff

            for protein in activities:
                for cmpd in activities[protein]:
                    for data in activities[protein][cmpd]:
                        act_value = data["act_value"]
                        pubmed_ids = data["pubmed_ids"]

                        if act_value >= cutoff:
                            # add protein
                            try:
                                _ = tmp[protein]
                            except KeyError:  # new protein
                                tmp[protein] = {}

                            # add cmpd for protein
                            try:
                                _ = tmp[protein][cmpd]
                            except KeyError:  # new cmpd for protein
                                tmp[protein][cmpd] = []

                            # add pubmed_ids and act_value for cmpd for protein
                            tmp[protein][cmpd].append(
                                {"act_value": act_value, "pubmed_ids": pubmed_ids}
                            )

        activities = tmp

    tmp = {}
    if by_year_threshold:
        for protein in activities:
            for cmpd in activities[protein]:
                for data in activities[protein][cmpd]:
                    act_value = data["act_value"]
                    pubmed_year = data["pubmed_ids"]

                    if pubmed_year is None:
                        continue

                    found = {}
                    for pubmed_id in pubmed_year:  # dict pubmed -> year
                        if pubmed_year[pubmed_id] is not None:
                            if pubmed_year[pubmed_id] <= by_year_threshold:
                                found[pubmed_id] = pubmed_year[pubmed_id]
                    if len(found):
                        # add protein
                        try:
                            _ = tmp[protein]
                        except KeyError:  # new protein
                            tmp[protein] = {}

                        # add cmpd for protein
                        try:
                            _ = tmp[protein][cmpd]
                        except KeyError:  # new cmpd for protein
                            tmp[protein][cmpd] = []

                        # add pubmed_ids and act_value for cmpd for protein
                        tmp[protein][cmpd].append(
                            {"act_value": act_value, "pubmed_ids": found}
                        )

        activities = tmp

    return activities


def getDTOforProtein(uniprot=None):
    """
    Get DTO id for protein.

    Args:
        uniprot (str/list/set/dict):    protein uniprot id

    Returns:
        dict:   {uniprot: dtoid}
    """

    results = {
        t: info["dtoid"] for t, info in getInfoForProtein(uniprot=uniprot).items()
    }
    return results


def getInfoForDTO(dtoid=None):
    """
    Get info for DTO id.

    Args:
        dtoid (str/list/set/dict):      dto id

    Returns:
        dict:   {dtoid: {'name': name, 'parent_id': parent_id, 'def': def, ...}}
    """

    dtoids = dtoid

    SELECT = ["dtoid", "name", "parent_id", "def"]
    FROM = ["dto"]

    if dtoids is not None:
        if isinstance(dtoids, str):
            dtoids = set([dtoids])
        dtoids = set([dtoid.replace("_", ":") for dtoid in dtoids if dtoid is not None])
        if not len(dtoids):  # all are None
            return {None: {"def": None, "parent_id": None, "name": None}}

        WHERE = ["dtoid IN ('{}')".format("','".join(list(dtoids)))]

    else:  # get info for all dtoids
        WHERE = ["1 = 1"]

    sql = "SELECT DISTINCT {} FROM {} WHERE {};".format(
        ",".join(SELECT), ",".join(FROM), " AND ".join(WHERE)
    )

    results = {}
    conn = openConnection()
    cur = conn.cursor(pymysql.cursors.DictCursor)
    cur.execute(sql)
    for row in cur:
        parent_id = row["parent_id"]
        if parent_id is not None:
            parent_id = parent_id.replace(":", "_")
        results[row["dtoid"]] = {
            "parent_id": parent_id,
            "name": row["name"],
            "def": row["def"],
        }

    cur.close()
    conn.close()

    for dtoid in dtoids.difference(set(results)):
        results[dtoid] = {"name": None, "parent_id": None, "def": None}

    results = {dtoid.replace(":", "_"): results[dtoid] for dtoid in results}

    return results


def getParentDTOforDTO(dtoid):
    dto_info = getInfoForDTO(dtoid=dtoid)
    dto_parent = {dtoid: dto_info[dtoid]["parent_id"] for dtoid in dto_info}

    return dto_parent


def getFullDTOforProtein(uniprot):
    if isinstance(uniprot, str):
        uniprot = set([uniprot])

    protein_dtoid = getDTOforProtein(uniprot=uniprot)

    dtoids = set([dtoid for dtoid in protein_dtoid.values() if dtoid is not None])
    while len(dtoids):
        parent_dtoid = getParentDTOforDTO(dtoid=dtoids)
        protein_dtoid = {**protein_dtoid, **parent_dtoid}
        dtoids = set([dtoid for dtoid in parent_dtoid.values() if dtoid is not None])

    results = {}
    for protein in uniprot:
        full_dtoid = []
        dtoid = protein_dtoid[protein]
        while dtoid is not None:
            full_dtoid.append(dtoid)
            dtoid = protein_dtoid[dtoid]

        if len(full_dtoid):
            full_dtoid = ".".join(full_dtoid[::-1])
        else:
            full_dtoid = None
        results[protein] = full_dtoid

    return results


def getInfoForPcid(pcid=None):
    """
    Get information for PANTHER protein class.

    Args:
        pcid (str/list/dict/set): protein panther class

    Returns:
        dict:   {pcid: {'name': name, 'parent_pcids': set of parent pcids, 'description': description}}

    """

    if isinstance(pcid, str):
        pcid = set([pcid])

    SELECT = ["pc.pcid", "pc.name", "pc.parent_pcids", "pc.description"]
    FROM = ["panther_class pc"]
    WHERE = ["1 = 1"]

    if pcid is not None:
        if isinstance(pcid, str):
            pcid = [pcid]
        WHERE.append("pc.pcid IN ('{}')".format("','".join(list(pcid))))

    sql = "SELECT DISTINCT {} FROM {} WHERE {};".format(
        ",".join(SELECT), ",".join(FROM), " AND ".join(WHERE)
    )

    results = {}
    conn = openConnection()
    cur = conn.cursor()
    cur.execute(sql)
    for pcid, name, parent_pcids, description in cur:
        if parent_pcids is None:
            parent_pcids = set()
        else:
            parent_pcids = set(parent_pcids.split("|"))
        results[pcid] = {
            "name": name,
            "parent_pcids": parent_pcids,
            "description": description,
        }

    cur.close()
    conn.close()

    return results


def getParentForPcid(pcid=None):
    pcid_info = getInfoForPcid(pcid=pcid)
    results = {pcid: pcid_info[pcid]["parent_pcids"] for pcid in pcid_info}
    return results


def getDeepestPcidForProtein(uniprot=None):
    """
    Get Panther class ID for protein.

    Args:
        uniprot (str/list/set/dict):    protein uniprot id

    Returns:
        dict:   {protein_id: {'pcid': pcid}}
    """

    # Get PCIDs for protein
    SELECT = ["p.uniprot", "pc.pcid"]
    FROM = ["protein p"]
    INNER_JOIN = [
        "p2pc ON p2pc.protein_id = p.id",
        "panther_class pc ON p2pc.panther_class_id = pc.id",
    ]
    WHERE = ["1 = 1"]

    if uniprot is not None:
        if isinstance(uniprot, str):
            uniprot = [uniprot]
        WHERE.append("p.uniprot IN ('{}')".format("','".join(list(uniprot))))

    sql = "SELECT DISTINCT {} FROM {} INNER JOIN {} WHERE {};".format(
        ",".join(SELECT),
        ",".join(FROM),
        " INNER JOIN ".join(INNER_JOIN),
        " AND ".join(WHERE),
    )

    results = {}
    conn = openConnection()
    cur = conn.cursor()
    cur.execute(sql)
    for uniprot, pcid in cur:
        try:
            results[uniprot].update([pcid])
        except KeyError:
            results[uniprot] = set([pcid])

    cur.close()
    conn.close()

    # save parent PCIDs relationships
    pcids = set([pcid for protein, pcids in results.items() for pcid in pcids])
    pcid_parents = getParentForPcid(pcid=pcids)
    parents = set([parent for parents in pcid_parents.values() for parent in parents])
    while len(parents):
        new_pcid_parents = getParentForPcid(pcid=parents)
        parents = set(
            [parent for parents in new_pcid_parents.values() for parent in parents]
        )
        pcid_parents = {**pcid_parents, **new_pcid_parents}  # update values

    # get full largest pcid hierarchical code for proteins
    results2 = {}
    for protein in results:
        pcids = results[protein]
        size = 0
        while len(pcids) > size:
            size = len(pcids)
            new_pcids = set()
            for pcid in pcids:
                outest_pcid = pcid.split(".")[
                    -1
                ]  # most superior/outer pcid in terms of hierarchy
                parent_pcids = pcid_parents[outest_pcid]
                for parent_pcid in parent_pcids:
                    new_pcids.update([pcid + "." + parent_pcid])
            pcids.update(new_pcids)

        largest_code = ""
        for pcid in pcids:
            if pcid.count(".") > largest_code.count("."):
                largest_code = pcid
        largest_code = ".".join(largest_code.split(".")[::-1])
        results2[protein] = largest_code

    return results2
