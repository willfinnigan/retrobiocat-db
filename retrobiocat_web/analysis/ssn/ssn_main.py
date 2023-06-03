from retrobiocat_web.analysis.ssn.ssn_save_load.ssn_save_load_sqlite import save_to_sqlite, load_edgelist_from_sqlite, \
    load_attrs_from_sqlite
from retrobiocat_web.mongo.models.biocatdb_models import Sequence, EnzymeType, UniRef50, SSN_record
from retrobiocat_web.analysis.sequence_analysis.all_by_all_blast import AllByAllBlaster
import mongoengine as db
import networkx as nx
import time
import json
import os
import pandas as pd
import shutil
from pathlib import Path
from retrobiocat_web.analysis import analysis_paths

class SSN(object):

    def __init__(self, enzyme_type, aba_blaster=None, log_level=0):

        self.graph = nx.Graph()

        self.min_score = 30

        self.enzyme_type = enzyme_type
        self.enzyme_type_obj = EnzymeType.objects(enzyme_type=enzyme_type)[0]

        if aba_blaster is None:
            self.aba_blaster = AllByAllBlaster(enzyme_type, log_level=log_level)
        else:
            self.aba_blaster = aba_blaster

        self.node_metadata = {}

        self.log_level = log_level

        self.ssn_folder = analysis_paths.ssn_folder

        self.save_path = f'{self.ssn_folder}/{self.enzyme_type}'
        Path(self.save_path).mkdir(parents=True, exist_ok=True)
        self.sqlite_file = f"{self.save_path}/{self.enzyme_type}_sqlite_ssn.sqlite"

        self.pos_path = f'{self.ssn_folder}/{self.enzyme_type}/positions'
        Path(self.pos_path).mkdir(parents=True, exist_ok=True)

        self.db_object = self._get_db_object()

        self.log(f"SSN object initialised for {enzyme_type}")

    def save(self):

        t0 = time.time()

        df_graph = nx.to_pandas_edgelist(self.graph)

        att_dict = {}
        for node in list(self.graph):
            att_dict[node] = self.graph.nodes[node]

        df_graph.to_csv(f"{self.save_path}/graph.csv", index=False)

        with open(f'{self.save_path}/attributes.json', 'wb') as outfile:
            outfile.write(json.dumps(att_dict).encode("utf-8"))

        t1 = time.time()

        self.log(f"Saved SSN as csv for {self.enzyme_type} in {round(t1 - t0, 1)} seconds")

    def save_sqlite(self):
        t0 = time.time()
        Path(self.save_path).mkdir(parents=True, exist_ok=True)
        save_to_sqlite(self.graph, self.sqlite_file)
        t1 = time.time()
        self.log(f"Saved SSN as sqlite for {self.enzyme_type} in {round(t1 - t0, 1)} seconds")

    def load_sqlite(self, alignment_score=0, node_list=None, include_no_edges=True, include_mutants=True, only_biocatdb=False):

        t0 = time.time()
        if not os.path.exists(self.sqlite_file):
            self.log(f"No saved SSN found for {self.enzyme_type}, could not load")
            return False

        edgelist = load_edgelist_from_sqlite(self.sqlite_file, alignment_score, node_list=node_list)
        attr_dict = load_attrs_from_sqlite(self.sqlite_file, node_list=node_list)
        self.graph = nx.from_pandas_edgelist(edgelist, edge_attr=['weight', 'i'])

        if include_no_edges == True:
            # Nodes with no edges are not in edge list..
            all_nodes = list(self.graph.nodes)
            for node in list(attr_dict.keys()):
                if node not in all_nodes:
                    self._add_protein_node(node)

        if include_mutants is False:
            self.filter_out_mutants()
        if only_biocatdb is True:
            self.filer_out_uniref()

        nx.set_node_attributes(self.graph, attr_dict)

        t3 = time.time()
        self.log(f"Loaded SSN from sqlite for {self.enzyme_type} in {round(t3 - t0, 1)} seconds")

    def load(self, include_mutants=True, only_biocatdb=False):

        t0 = time.time()
        if not os.path.exists(f"{self.save_path}/graph.csv") or not os.path.exists(f"{self.save_path}/attributes.json"):
            self.log(f"No saved SSN found for {self.enzyme_type}, could not load")
            return False

        #df_graph = pd.read_csv(f"{self.save_path}/graph.csv")

        df_graph = pd.read_csv(f"{self.save_path}/graph.csv")
        att_dict = json.load(open(f'{self.save_path}/attributes.json'))
        t1 = time.time()

        self.graph = nx.from_pandas_edgelist(df_graph, edge_attr=['weight', 'i'])

        t2 = time.time()

        # Nodes with no edges are not in edge list..
        for node in list(att_dict.keys()):
            if node not in list(self.graph.nodes):
                self._add_protein_node(node)

        if include_mutants is False:
            self.filter_out_mutants()
        if only_biocatdb is True:
            self.filer_out_uniref()

        nx.set_node_attributes(self.graph, att_dict)

        t3 = time.time()
        self.log(f"Loaded SSN from csv for {self.enzyme_type} in {round(t3 - t0, 1)} seconds")
        self.log(f"- {round(t1 - t0, 1)} seconds to load dataframes")
        self.log(f"- {round(t2 - t1, 1)} seconds to load actual ssn from edge list")
        self.log(f"- {round(t3 - t2, 1)} seconds for attributes")

    def add_protein(self, seq_obj):
        """ Add the protein to the graph, along with any proteins which have alignments """

        self.log(f"Adding node - {seq_obj.enzyme_name} and making alignments..")
        t0 = time.time()

        name = seq_obj.enzyme_name
        self._add_protein_node(name, alignments_made=True)
        alignment_names, alignment_scores, identities, coverages = self.aba_blaster.get_alignments(seq_obj)
        # self.graph.nodes[name]['attributes']['alignments_made'] = True

        count = 0
        for i, protein_name in enumerate(alignment_names):
            count += self._add_protein_node(protein_name)
            self._add_alignment_edge(seq_obj.enzyme_name, protein_name, alignment_scores[i], identities[i],
                                     coverages[i])

        t1 = time.time()
        seq_obj.alignments_made = True
        seq_obj.save()

        self.log(f"{count} new nodes made for alignments, with {len(alignment_names)} edges added")
        self.log(f"Protein {seq_obj.enzyme_name} processed in {round(t1 - t0, 0)} seconds")

    def add_multiple_proteins(self, list_seq_obj):
        for seq_obj in list_seq_obj:
            self.add_protein(seq_obj)

    def nodes_need_alignments(self, max_num=None):
        """ Return nodes which needs alignments making, maximum of max_num"""

        t0 = time.time()
        need_alignments = []
        count = 0
        for node in list(self.graph.nodes):
            if dict(self.graph.nodes[node]).get('alignments_made', False) == False:
                seq_obj = self._get_sequence_object(node)
                need_alignments.append(seq_obj)
                count += 1
                if count == max_num:
                    break

        t1 = time.time()
        self.log(f"Identified {count} nodes which need alignments making in {round(t1 - t0, 1)} seconds")

        return need_alignments

    def nodes_not_present(self, only_biocatdb=False, max_num=None):
        """ Return a list of enzymes which are not in the ssn """

        # Get a list of all sequence objects of enzyme type
        t0 = time.time()
        sequences = Sequence.objects(db.Q(enzyme_type=self.enzyme_type) &
                                     db.Q(sequence__ne="") &
                                     db.Q(sequence__ne=None) &
                                     db.Q(sequence_unavailable__ne=True) &
                                     db.Q(reviewed=True))
        if only_biocatdb is True:
            seq_objects = list(sequences)
        else:
            unirefs = UniRef50.objects(enzyme_type=self.enzyme_type_obj)
            seq_objects = list(sequences) + list(unirefs)

        # Get sequences not in nodes
        not_in_nodes = []
        for seq_obj in seq_objects:
            if seq_obj.enzyme_name not in list(self.graph.nodes):
                if seq_obj.sequence != None:
                    if len(seq_obj.sequence) > 12:
                        not_in_nodes.append(seq_obj)

        # Return only up to the maximum number of sequences
        if max_num != None:
            if len(not_in_nodes) > max_num:
                not_in_nodes = not_in_nodes[0:max_num]

        t1 = time.time()
        self.log(
            f"Identified {len(not_in_nodes)} {self.enzyme_type} proteins which need adding, in {round(t1 - t0, 1)} seconds")
        return not_in_nodes

    def remove_nonexisting_seqs(self):

        t0 = time.time()
        sequences = Sequence.objects(enzyme_type=self.enzyme_type).distinct('enzyme_name')
        unirefs = UniRef50.objects(enzyme_type=self.enzyme_type_obj).distinct('enzyme_name')
        protein_names = list(sequences) + list(unirefs)
        count = 0
        for node in list(self.graph.nodes):
            if node not in protein_names:
                self.log(f"Node: {node} not in the database - removing")
                self.graph.remove_node(node)
                count += 1

        t1 = time.time()
        self.log(f"Identified {count} sequences which were in SSN but not in database, in {round(t1 - t0, 1)} seconds")

    def remove_seqs_marked_with_no_alignments(self):

        t0 = time.time()
        enz_type_q = db.Q(enzyme_type=self.enzyme_type)
        align_q = db.Q(alignments_made__ne=True)
        protein_names = Sequence.objects(enz_type_q & align_q).distinct('enzyme_name')

        count = 0
        for node in list(self.graph.nodes):
            if node in protein_names:
                self.log(f"Node: {node} marked as having no alignments made - removing")
                self.graph.remove_node(node)
                count += 1

        t1 = time.time()
        self.log(f"Identified {count} sequences which were in SSN but not are marked has not having their alignments made, in {round(t1 - t0, 1)} seconds")

    def get_graph_filtered_edges(self, alignment_score, min_ident=0):

        sub_graph = nx.Graph([(u, v, d) for u, v, d in self.graph.edges(data=True) if d['weight'] >= alignment_score and d['i'] >= min_ident])
        for node in self.graph.nodes:
            if node not in sub_graph.nodes:
                sub_graph.add_node(node)

        return sub_graph

    def delete_edges_below_alignment_score(self, alignment_score, min_ident=0.0):
        self.log(f"Removing edges below alignment score: {alignment_score}, min identity: {min_ident}")
        edge_list_to_remove = []
        for edge in self.graph.edges(data=True):
            if edge[2]['weight'] < alignment_score or edge[2]['weight'] < min_ident:
                edge_list_to_remove.append(edge)

        for edge in edge_list_to_remove:
            self.graph.remove_edge(edge[0], edge[1])

    def filter_out_mutants(self):
        t0 = time.time()
        mutants = Sequence.objects(db.Q(enzyme_type=self.enzyme_type) &
                                   (db.Q(mutant_of__ne='') & db.Q(mutant_of__ne=None))).distinct('enzyme_name')

        for mutant in list(mutants):
            if mutant in self.graph.nodes:
                self.graph.remove_node(mutant)

        t1 = time.time()
        self.log(f'Filtered mutants from graph in {round(t1 - t0, 1)} seconds')

    def filer_out_uniref(self):
        t0 = time.time()
        for node in list(self.graph.nodes):
            if 'UniRef50' in node:
                self.graph.remove_node(node)

        t1 = time.time()
        self.log(f'Filtered uniref50 sequences from graph in {round(t1 - t0, 1)} seconds')

    def only_specificed_nodes(self, spec_nodes):
        t0 = time.time()
        for node in list(self.graph.nodes):
            if node not in spec_nodes:
                self.graph.remove_node(node)
        t1 = time.time()
        self.log(f'Filtered to only specified sequences in {round(t1 - t0, 1)} seconds')


    def _add_protein_node(self, node_name, alignments_made=False):
        """ If a protein is not already in the graph, then add it """
        if 'UniRef50' in node_name:
            node_type = 'uniref'
        else:
            node_type = 'biocatdb'

        if node_name not in self.graph.nodes:
            self.graph.add_node(node_name, node_type=node_type,
                                alignments_made=alignments_made)
            return 1

        if alignments_made == True:
            self.graph.nodes[node_name]['alignments_made'] = True

        return 0

    def _add_alignment_edge(self, node_name, alignment_node_name, alignment_score, i, c):
        if node_name != alignment_node_name:
            if alignment_score > self.min_score:
                self.graph.add_edge(node_name, alignment_node_name, weight=alignment_score, i=i)

    def log(self, msg, level=10):
        if level >= self.log_level:
            print("SSN: " + msg)

    @staticmethod
    def _get_sequence_object(enzyme_name):
        if 'UniRef50' in enzyme_name:
            return UniRef50.objects(enzyme_name=enzyme_name)[0]
        else:
            return Sequence.objects(enzyme_name=enzyme_name)[0]

    def _get_db_object(self):
        """ Either finds existing db entry for ssn of enzyme type, or makes a new one """

        query = SSN_record.objects(enzyme_type=self.enzyme_type_obj)
        self.log(f"{len(list(query))} ssn objects for {self.enzyme_type}", level=3)
        if len(query) == 0:
            db_ssn = SSN_record(enzyme_type=self.enzyme_type_obj)
        else:
            db_ssn = query[0]

        return db_ssn

    def set_status(self, status):
        self.db_object.status = status
        self.db_object.save()

    def get_clusters(self, alignment_score):
        t0 = time.time()
        graph = self.get_graph_filtered_edges(alignment_score)
        clusters = list(nx.connected_components(graph))
        with_uniref = []
        without_uniref = []

        for cluster in clusters:
            new_cluster_with = []
            new_cluster_without = []
            for node in cluster:
                new_cluster_with.append(node)
                if 'UniRef50_' not in node:
                    new_cluster_without.append(node)

            if len(new_cluster_with) != 0:
                with_uniref.append(new_cluster_with)
            if len(new_cluster_without) != 0:
                without_uniref.append(new_cluster_without)

        t1 = time.time()
        self.log(f"Retrieved clusters for {self.enzyme_type} at alignment score {alignment_score} in {round(t1 - t0, 1)} seconds")

        return without_uniref, with_uniref

    def clear_position_information(self):
        self.log(f"Clearing old node position information", level=1)
        if self.db_object.precalc_status != 'Empty':
            self.db_object.num_at_alignment_score = {}
            self.db_object.pos_at_alignment_score = {}
            self.db_object.precalculated_vis = {}
            self.db_object.identity_at_alignment_score = {}
            self.db_object.precalc_status = 'Empty'
            self.db_object.save()

            shutil.rmtree(self.pos_path)
            Path(self.pos_path).mkdir(parents=True, exist_ok=True)

    def clear_all_data(self):
        self.clear_position_information()
        shutil.rmtree(self.save_path, ignore_errors=True)
        shutil.rmtree(self.pos_path, ignore_errors=True)
        self.db_object.delete()

    def rename_dir(self, new_folder_name):
        rename_from = self.save_path
        rename_to = f'{self.ssn_folder}/{new_folder_name}'
        shutil.move(rename_from, rename_to, copy_function=shutil.copytree)
        print(f"renamed ssn folder to {new_folder_name}")

    def send_to_cytoscape(self, alignment_score, selected_nodes=None):

        # filtering
        graph = self.get_graph_filtered_edges(alignment_score, min_ident=0)
        if selected_nodes is not None:
            for node in selected_nodes:
                if node in graph.nodes:
                    graph.remove_node(node)

        # add uniref metadata
        unirefs = UniRef50.objects(enzyme_type=self.enzyme_type_obj).exclude('id', 'enzyme_type', 'sequence',
                                                                             "result_of_blasts_for", 'alignments_made')
        for seq_obj in unirefs:
            update_dict = json.loads(seq_obj.to_json())
            if seq_obj.enzyme_name in graph.nodes:
                graph.nodes[seq_obj.enzyme_name]['pfam_codes'] = update_dict.get('pfam_codes', [])
                graph.nodes[seq_obj.enzyme_name]['pfams'] = update_dict.get('pfams', {})
                graph.nodes[seq_obj.enzyme_name]['sp_annotated'] = update_dict.get('sp_annotated', None)
                graph.nodes[seq_obj.enzyme_name]['num_uni90'] = update_dict.get('num_uni90', 0)
                graph.nodes[seq_obj.enzyme_name]['num_uni100'] = update_dict.get('num_uni100', 0)
                graph.nodes[seq_obj.enzyme_name]['rep_name'] = update_dict.get('rep_name', '')
                graph.nodes[seq_obj.enzyme_name]['protein_name'] = update_dict.get('protein_name', '')
                graph.nodes[seq_obj.enzyme_name]['tax'] = update_dict.get('tax', '')
                graph.nodes[seq_obj.enzyme_name]['tax_id'] = update_dict.get('tax_id', '')
                graph.nodes[seq_obj.enzyme_name]['rep_pdbs'] = update_dict.get('rep_pdbs', [])

        # convert graph to cytoscape format
        cytoscape = nx.cytoscape_data(graph)

        return cytoscape

if __name__ == "__main__":
    print(str(Path(__file__).parents[1]) + f'/analysis_data/ssn')
    from retrobiocat_web.mongo.default_connection import make_default_connection
    make_default_connection()

    from retrobiocat_web.analysis.ssn.ssn_vis import SSN_Visualiser



    enzyme_type = 'TA'

    aba_blaster = AllByAllBlaster(enzyme_type, log_level=1)
    aba_blaster.make_blast_db()
    ssn = SSN(enzyme_type, aba_blaster=aba_blaster)
    ssn.load_sqlite()

    biocatdb_seqs = ssn.nodes_not_present(only_biocatdb=True, max_num=200)
    print(biocatdb_seqs)

    ssn.add_multiple_proteins(biocatdb_seqs)

    biocatdb_seqs = ssn.nodes_not_present(only_biocatdb=True, max_num=200)
    print(biocatdb_seqs)

    df_graph = nx.to_pandas_edgelist(ssn.graph)
    print(df_graph.tail())
    print(df_graph[df_graph['target'] == 'BmTA_S119G'])
    print(df_graph[df_graph['target'] == 'ArRmut11- oTA'])
    print(df_graph[df_graph['target'] == 'VfTA_L56A-I259T'])

    ssn.save_sqlite()

    all_nodes = list(ssn.graph.nodes)


    #ssn.graph.add_edge('test_node', 'ARmutTA')

    #ssn.save_sqlite()

    #nodes = [n for n in ssn.graph.nodes if ssn.graph.nodes[n]['node_type'] == 'biocatdb']
    #print(nodes)

    #att_dict = {}
    #for node in list(ssn.graph):
    #    att_dict[node] = ssn.graph.nodes[node]
    #print(att_dict)
    #ssn.save_sqlite()

    #vis = SSN_Visualiser(enzyme_type, log_level=1)
    #nodes, edges = vis.visualise(ssn, 46)
    #cytoscape = ssn.send_to_cytoscape(45)
    #with open('data.cyjs', 'w') as outfile:
    #    json.dump(cytoscape, outfile)

