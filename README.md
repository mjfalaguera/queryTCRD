# queryTCRD

queryTCRD is an (unofficial) python module with functions to query the Target Central Resource Database (TCRD).

### Usage
1) Download and install the latest SQL version of the TCRD from http://juniper.health.unm.edu/tcrd/download/.
2) Edit your database details in the first rows of the queryTCRD.py script.
3) Install pymysql package.
3) Import queryTCRD.py as a module in your scripts and start playing!

### Extras
Find a detailed view of the TCRD schema in TCRDv6_ERD.pdf (source: http://juniper.health.unm.edu/tcrd/download/TCRDv6_ERD.pdf) and suggest new functions to add.

### What is the TCRD?

(Adapted from http://juniper.health.unm.edu/tcrd/)

The **Target Central Resource Database (TCRD)** is the central resource behind the Illuminating the Druggable Genome Knowledge Management Center (IDG-KMC). TCRD contains information about human targets, with special emphasis on four families of targets that are central to the NIH IDG initiative: GPCRs, kinases and ion channels (note that olfactory GPCRs are treated as a separate family. A key aim of the KMC is to classify the development/druggability level of targets. We currently categorize targets into four **target development/druggability levels (TDLs)** defined as follows:
- **Tclin**, these targets have activities in DrugCentral (ie. approved drugs) with known mechanism of action.
- **Tchem**, these targets have activities in ChEMBL, Guide to Pharmacology or DrugCentral that satisfy IDG activity thresholds.
- **Tbio**, non-Tchem/Tclin targets that satisfy one or more of the following criteria: target is above the cutoff criteria for Tdark and is annotated with a Gene Ontology Molecular Function or Biological Process leaf term(s) with an Experimental Evidence code.
- **Tdark**, these are targets about which virtually nothing is known and satisfy two or more of the following criteria: a PubMed text-mining score from Jensen Lab < 5, <= 3 Gene RIFs, <= 50 Antibodies available according to antibodypedia.com.
