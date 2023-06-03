import pytest
from retrobiocat_web.analysis.uniprot.uniref_data_extraction import check_id_match, get_uniref_members
from retrobiocat_web.analysis.uniprot.xml_retrieval import UniRef_XML_Retriever

uniref_codes_to_test = ["UniRef50_Q9C969"]

@pytest.fixture
def xml_dict_code(request):
    code = request.param
    uniref_xml_retriever = UniRef_XML_Retriever(log_level='DEBUG')
    xml_dict = uniref_xml_retriever.get_xml(code)
    return [xml_dict, code]

class Test_UniRef_Data_from_XML():

    @pytest.mark.parametrize('xml_dict_code', uniref_codes_to_test, indirect=True)
    def test_get_cluster_id_matches(self, xml_dict_code):
        xml_dict = xml_dict_code[0]
        code = xml_dict_code[1]
        assert check_id_match(code, xml_dict) is True

    @pytest.mark.parametrize('xml_dict_code', uniref_codes_to_test, indirect=True)
    def test_has_uniref_90_members(self, xml_dict_code):
        xml_dict = xml_dict_code[0]
        uni90, uni100, uniprot, uniprot_kb = get_uniref_members(xml_dict)
        assert len(uni90) > 0

    @pytest.mark.parametrize('xml_dict_code', uniref_codes_to_test, indirect=True)
    def test_has_uniref_100_members(self, xml_dict_code):
        xml_dict = xml_dict_code[0]
        uni90, uni100, uniprot, uniprot_kb = get_uniref_members(xml_dict)
        assert len(uni100) > 0

    @pytest.mark.parametrize('xml_dict_code', uniref_codes_to_test, indirect=True)
    def test_has_uniprot_members(self, xml_dict_code):
        xml_dict = xml_dict_code[0]
        uni90, uni100, uniprot, uniprot_kb = get_uniref_members(xml_dict)
        assert len(uniprot) > 0



