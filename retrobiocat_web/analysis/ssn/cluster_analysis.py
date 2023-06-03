from retrobiocat_web.analysis.ssn import ssn_main

'''
First make SSN, then get clusters.
Then do an analysis of the clusters
'''

def get_ssn(enzyme_type, alignment_score, ):
    ssn = ssn_main.SSN(enzyme_type)
    ssn.load_sqlite()

    clusters, clusters_no_uniref = ssn.get_clusters(alignment_score)

