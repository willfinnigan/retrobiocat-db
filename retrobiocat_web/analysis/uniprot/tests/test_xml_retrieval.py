from collections import OrderedDict
import pytest
from retrobiocat_web.analysis.uniprot.xml_retrieval import XML_Retriever, UniRef_XML_Retriever, UniProt_XML_Retriever

uniref_codes_to_test = ["UniRef50_Q9C969"]
uniprot_codes_to_test = ["Q9C969"]

@pytest.mark.parametrize('uniref_code', uniref_codes_to_test)
class Test_XML_UniRef_Retrieval():

    def test_can_get_an_uniref_xmls(self, uniref_code):
        uniref_retriever = UniRef_XML_Retriever(log_level='DEBUG')
        xml = uniref_retriever._retrieve_xml(uniref_code)
        assert xml is not None

    def test_can_parse_uniref_xml(self, uniref_code):
        uniref_retriever = UniRef_XML_Retriever(log_level='DEBUG')
        xml_dict = uniref_retriever.get_xml(uniref_code)
        assert isinstance(xml_dict, OrderedDict)

@pytest.mark.parametrize('uniprot_code', uniprot_codes_to_test)
class Test_XML_UniProt_Retrieval():

    def test_can_get_an_uniprot_xmls(self, uniprot_code):
        uniprot_retriever = UniProt_XML_Retriever(log_level='DEBUG')
        xml = uniprot_retriever._retrieve_xml(uniprot_code)
        print(xml)
        assert xml is not None

    def test_can_parse_uniprot_xml(self, uniprot_code):
        uniprot_retriever = UniProt_XML_Retriever(log_level='DEBUG')
        xml_dict = uniprot_retriever.get_xml(uniprot_code)
        assert isinstance(xml_dict, OrderedDict)


