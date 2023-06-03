from retrobiocat_web.mongo.functions.save_activity_data import cascade_activity_data, check_activity_data
from retrobiocat_web.mongo.model_queries import sequence_queries
from rdkit import Chem

class ActivityDataProcessor(object):
    """
    Class for processing a list of activity data so its ready for saving.
    Checks all the data is ok, and cascades activity data down the various levels
    eg. Conversion -> Categorical -> Binary
    """

    def __init__(self):
        """Default values for cascading activity data"""
        self.CATS_CONVERSION = {'High': 60,
                                'Medium': 10,
                                'Low': 0.1,
                                'None': 0}

        self.CATS_SA = {'High': 1,
                        'Medium': 0.25,
                        'Low': 0.01,
                        'None': 0}

        self.DEFAULT_CONC = 10

        self.smi_cols = ['substrate_1_smiles', 'substrate_2_smiles', 'product_1_smiles']

    @staticmethod
    def _add_row_numbers(data):
        """Adds row numbers, which are used by tabulator in web form"""
        for i, data_dict in enumerate(data):
            data[i]['n'] = i+1
        return data

    def _run_activity_cascade(self, data_dict):
        """Cascades activity data, for example getting categorical from conversion ect.."""

        data_dict = cascade_activity_data.sa_from_kinetics(data_dict, self.DEFAULT_CONC)
        data_dict = cascade_activity_data.category_from_sa(data_dict, self.CATS_SA)
        data_dict = cascade_activity_data.category_from_conversion(data_dict, self.CATS_CONVERSION)
        data_dict = cascade_activity_data.binary_from_category(data_dict)
        data_dict = cascade_activity_data.make_sure_binary_is_int(data_dict)
        return data_dict

    @staticmethod
    def _check_data_dict(data_dict, paper):
        """Checks that data conforms to requirements for activity data"""
        issues = []
        issues += check_activity_data.check_required_columns(data_dict)
        issues += check_activity_data.check_all_required_are_numbers(data_dict)
        issues += check_activity_data.check_all_have_binary(data_dict)
        paper_seqs = sequence_queries.seqs_of_paper(paper)
        issues += check_activity_data.check_seqs_are_defined(data_dict, paper_seqs)
        return issues

    def _convert_to_rdkit_smiles(self, data_dict):
        """Ensures that every smiles is the default RDKit smiles"""

        issues = []
        for col in self.smi_cols:
            if col in data_dict:
                if data_dict[col] != '' or data_dict[col] is not None:
                    smi = data_dict[col]
                    try:
                        mol = Chem.MolFromSmiles(smi)
                        smi = Chem.MolToSmiles(mol)
                    except:
                        issues.append(f"Row {data_dict['n']}, column {col} - smi {smi} was not accepted by RdKit")
                    data_dict[col] = smi

        return issues, data_dict

    def process_data(self, data, paper, add_row_numbers=False):
        """ Process an activity data list ready for saving """

        all_issues = []

        # remove empty
        data = [x for x in data if x != {}]

        if add_row_numbers:
            data = self._add_row_numbers(data)   # shouldn't need to add row numbers, they are added by tabulator

        # process the data
        processed_data = []
        for data_dict in data:
            smi_issues, data_dict = self._convert_to_rdkit_smiles(data_dict)  # ensure smis are accepted by rdkit
            data_dict = self._run_activity_cascade(data_dict)
            other_issues = self._check_data_dict(data_dict, paper)

            all_issues += other_issues
            all_issues += smi_issues


            if len(all_issues) == 0:
                processed_data.append(data_dict)

        return processed_data, all_issues
