import pytest

from retrobiocat_web.analysis.uniprot.uniprot_data_extraction import get_pfams, get_rhea, is_swissprot, get_pdbs
from retrobiocat_web.analysis.uniprot.xml_retrieval import UniProt_XML_Retriever

uniprot_codes_to_test = [
    ["Q9C969", {'expected_pfam': {'PF00155': 'Aminotran_1_2'},
                'expected_rhea': ['RHEA:31763', 'RHEA:14093', 'RHEA:15093'],
                'expected_pdbs': [],
                'expected_is_swiss': True},
     ]
]

@pytest.fixture
def xml_dict_code(request):
    code = request.param
    uniprot_xml_retriever = UniProt_XML_Retriever(log_level='DEBUG')
    xml_dict = uniprot_xml_retriever.get_xml(code)
    return [xml_dict, code]

@pytest.mark.parametrize('xml_dict_code, expected', uniprot_codes_to_test, indirect=['xml_dict_code'])
class Test_UniProt_Data_from_XML():

    def test_get_pfams(self, xml_dict_code, expected):
        xml_dict = xml_dict_code[0]
        pfams = get_pfams(xml_dict)
        print(pfams)
        assert pfams == expected['expected_pfam']

    def test_get_rhea(self, xml_dict_code, expected):
        xml_dict = xml_dict_code[0]
        rhea = get_rhea(xml_dict)
        print(rhea)
        assert rhea == expected['expected_rhea']

    def test_get_pdbs(self, xml_dict_code, expected):
        xml_dict = xml_dict_code[0]
        pdbs = get_pdbs(xml_dict)
        assert pdbs == expected['expected_pdbs']

    def test_is_swiss_prot(self, xml_dict_code, expected):
        xml_dict = xml_dict_code[0]
        is_swiss = is_swissprot(xml_dict)
        print(is_swiss)
        assert is_swiss == expected['expected_is_swiss']

