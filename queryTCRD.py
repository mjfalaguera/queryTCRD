#!/usr/bin/env python

"""
queryTCRD.py: Queryng the Target Central Resource Database (TCRD).
"""

__author__ = "Maria J. Falaguera"
__date__ = "22 Jun 2023"

import pymysql

# Set your DB details
dbname = "tcrdv6_12_4"
user = ""
host = ""
password = ""


def openConnection(dbname=dbname):
    """
    Open database connection
    """
    schema = "public"
    conn = pymysql.connect(
        database=dbname,
        user=user,
        host=host,
        password=password,
    )
    return conn


def getProteins():
    """
    Get all proteins (from protein table).
    """

    SELECT = ["protein.id", "protein.uniprot"]
    FROM = ["protein"]

    sql = "SELECT DISTINCT {} FROM {};".format(",".join(SELECT), ",".join(FROM))

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
    """
    Count proteins.
    """
    return len(getProteins())


def getTargets():
    """
    Get all targets (from target table).
    """

    SELECT = ["target.id"]
    FROM = ["target"]

    sql = "SELECT DISTINCT {} FROM {};".format(
        ",".join(SELECT),
        ",".join(FROM),
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


def getInfoForProtein(uniprot=None):
    """
    Get info for protein.

    Args:
        uniprot (str/list/set/dict):    protein uniprot id
        idx (str/list/set/dict):        protein id in TCRD

    Returns:
        dict:   {uniprot: {'tdl': tdl, 'name': name, 'family': family, ...}}
    """

    SELECT = [
        "p.uniprot",
        "p.name",
        "t.tdl",
        "t.fam as target_family",
        "p.sym as gene_name",
        "p.family as protein_family",
        "t.ttype",
        "p.dtoid",
        "p.description",
        "p.dtoclass",
        "p.id as protein_id",
        "t.id as target_id",
    ]
    FROM = ["protein p"]
    INNER_JOIN = [
        "t2tc ON t2tc.protein_id = p.id",
        "target t ON t2tc.target_id = t.id",
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
    cur = conn.cursor(pymysql.cursors.DictCursor)
    cur.execute(sql)
    for row in cur:
        results[row["uniprot"]] = {r: row[r] for r in row if r != "uniprot"}

    for uniprot in results:
        if results[uniprot]["protein_family"] is None:
            results[uniprot]["protein_family"] = ""
        if results[uniprot]["target_family"] is None:
            results[uniprot]["target_family"] = ""

    cur.close()
    conn.close()

    return results


def getTDLforProtein(uniprot=None):
    """
    Get TDL for protein.

    Args:
        uniprot (str/list/set/dict):    protein uniprot id

    Returns:
        dict:   {uniprot: TDL}
    """

    results = {t: info["tdl"] for t, info in getInfoForProtein(uniprot=uniprot).items()}
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


def getChemblRefYear(uniprot=None, include_empty=False):
    """
    Get ChEMBL First Reference Year for protein.

    Args:
        uniprot (str/list/set/dict):    protein uniprot id
        include_empty (bool):           include proteints with no ChEMBL First Reference Year annotated

    Returns:
        dict:   {uniprot: chembl_ref_year, ...}}
    """

    SELECT = ["p.uniprot", "tdl_info.integer_value as chembl_ref_year"]
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
        results[row["uniprot"]] = row["chembl_ref_year"]

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
    Get IDG target family for protein.

    Args:
        uniprot (str/list/set/dict):    protein uniprot id

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


def getCmpdActivitiesForProtein(uniprot=None, smiles=False, pact_cutoff=False):
    """
    Get compound activities for protein.

    Args:
        uniprot (str/list/set/tuple/dict):    protein to query
        smiles (bool):                        only activities implying cmpd with SMILES
        pact_cutoff (bool/float/int):        only activities with a value above cutoff

    Returns:
        dict: { uniprot: {'cmpd_id': [{'act_value': act_value, 'pubmed_ids': pubmed_ids set}, ...], ... }
    """

    WHERE = ["1 = 1"]

    if smiles:
        WHERE.append("cmpd_act.smiles != ''")

    if pact_cutoff:
        WHERE.append("cmpd_act.act_value >= {}".format(pact_cutoff))

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
        "p.uniprot",
        "cmpd_act.cmpd_id_in_src",
        "cmpd_act.act_value",
        "cmpd_act.pubmed_ids",
    ]
    FROM = ["protein p"]
    INNER_JOIN = [
        "t2tc ON t2tc.protein_id = p.id",
        "cmpd_activity cmpd_act ON cmpd_act.target_id = t2tc.target_id",
    ]
    sql = "SELECT {} FROM {} INNER JOIN {} WHERE {};".format(
        ",".join(SELECT),
        ",".join(FROM),
        " INNER JOIN ".join(INNER_JOIN),
        " AND ".join(WHERE),
    )

    conn = openConnection()
    cur = conn.cursor()
    cur.execute(sql)

    cmpd_acts = {}
    for r in cur:
        protein, cmpd, act, pubmed_ids = r

        # add protein
        try:
            _ = cmpd_acts[protein]
        except KeyError:  # new protein
            cmpd_acts[protein] = {}

        # add cmpd for protein
        try:
            _ = cmpd_acts[protein][cmpd]
        except KeyError:  # new cmpd for protein
            cmpd_acts[protein][cmpd] = []

        # add pubmed_ids and act_value for cmpd for protein
        if pubmed_ids is not None:
            pubmed_ids = set(
                [int(pubmed_id) for pubmed_id in set(pubmed_ids.split("|"))]
            )

        cmpd_acts[protein][cmpd].append(
            {"act_value": float(act), "pubmed_ids": pubmed_ids}
        )

    cur.close()
    conn.close()

    # annotate pubmed_ids with their publication year
    tmp = set()
    for p in cmpd_acts:
        for cmpd in cmpd_acts[p]:
            for act in cmpd_acts[p][cmpd]:
                pubmed_ids = act["pubmed_ids"]
                if pubmed_ids is not None:
                    tmp.update(pubmed_ids)

    pubmed_year = getYearForPubmed(pubmed_id=tmp)

    tmp = {}
    for p in cmpd_acts:
        tmp[p] = {}
        for cmpd in cmpd_acts[p]:
            tmp[p][cmpd] = []
            for act in cmpd_acts[p][cmpd]:
                # assign year to pubmed_ids
                if act["pubmed_ids"] is None:
                    tmp[p][cmpd].append(
                        {"act_value": act["act_value"], "pubmed_ids": {None: None}}
                    )

                else:
                    pubmed_ids = {}
                    for pubmed_id in act["pubmed_ids"]:
                        try:
                            year = pubmed_year[pubmed_id]
                        except KeyError:  # no year found for pubmed_id
                            year = None
                        pubmed_ids[pubmed_id] = year
                    tmp[p][cmpd].append(
                        {"act_value": act["act_value"], "pubmed_ids": pubmed_ids}
                    )
    cmpd_acts = tmp

    return cmpd_acts


def getDrugActivitiesForProtein(uniprot=None, smiles=False, pact_cutoff=False):
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
        WHERE.append("drug_act.smiles IS NOT NULL")

    if pact_cutoff:
        WHERE.append("drug_act.act_value >= {}".format(pact_cutoff))

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

    # Get drug activities
    SELECT = ["p.uniprot", "drug_act.drug", "drug_act.act_value"]
    FROM = ["protein p"]
    INNER_JOIN = [
        "t2tc ON t2tc.protein_id = p.id",
        "drug_activity drug_act ON drug_act.target_id = t2tc.target_id",
    ]
    sql = "SELECT {} FROM {} INNER JOIN {} WHERE {};".format(
        ",".join(SELECT),
        ",".join(FROM),
        " INNER JOIN ".join(INNER_JOIN),
        " AND ".join(WHERE),
    )

    conn = openConnection()
    cur = conn.cursor()
    cur.execute(sql)

    drug_acts = {}
    for r in cur:
        protein, drug, act = r
        if act is not None:
            act = float(act)

        # add protein
        try:
            _ = drug_acts[protein]
        except KeyError:  # new protein
            drug_acts[protein] = {}

        # add cmpd for protein
        try:
            _ = drug_acts[protein][drug]
        except KeyError:  # new cmpd for protein
            drug_acts[protein][drug] = []

        drug_acts[protein][drug].append({"act_value": act, "pubmed_ids": set()})

    cur.close()

    return drug_acts


def getCutoffForFamily(family=None):
    """
    Dictionary of IDG family-specific pAct cutoffs.

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

    SELECT = ["cmpd_id_in_src", "smiles"]
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
    cur = conn.cursor()
    cur.execute(sql)
    for cmpd_id_in_src, smiles in cur:
        results[cmpd_id_in_src] = {"smiles": smiles}

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

    results = {c: info["smiles"] for c, info in getInfoForCmpd(cmpd=cmpd).items()}
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
    cur = conn.cursor()
    cur.execute(sql)
    for drug, smiles in cur:
        results[drug] = {"smiles": smiles}

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

    results = {d: info["smiles"] for d, info in getInfoForDrug(drug=drug).items()}
    return results
