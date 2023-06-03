import requests
import xmltodict
import time
from retrobiocat_web.logging import add_logger


class XML_Retriever():

    def __init__(self, base_url, retry_num=3, retry_wait=1, log_level='WARNING'):
        self.base_url = base_url
        self.retry_num = retry_num
        self.retry_wait = retry_wait
        self.logger = add_logger('XML_Retriever', level=log_level)

    def get_xml(self, identifier):
        xml = self._retrieve_xml(identifier)
        xml_dict = self._xml_to_dict(xml, identifier)
        return xml_dict

    def _xml_to_dict(self, xml, identifier):
        if xml is None:
            return {}

        try:
            xml_dict = xmltodict.parse(xml)
            self.logger.debug(f"XML dict parsed for {identifier}")
            return xml_dict

        except Exception as e:
            self.logger.debug(f"Exception getting xml for {identifier}")
            self.logger.debug(str(e))
            return {}

    def _retrieve_xml(self, identifier):
        url = f"{self.base_url}/{identifier}.xml"

        for i in range(self.retry_num):
            try:
                req = requests.get(url)

                if req.status_code in [200]:
                    xml = req.text
                    return xml
                else:
                    self.logger.debug(f"Failed to retrieve xml for {identifier}. Status code = {req.status_code}")
                    return None

            except:
                time.sleep(self.retry_wait)

        return None


class UniRef_XML_Retriever(XML_Retriever):

    def __init__(self, retry_num=3, retry_wait=1, log_level='WARNING'):
        self.base_url = "https://www.uniprot.org/uniref"
        self.retry_num = retry_num
        self.retry_wait = retry_wait
        self.logger = add_logger('UniRef_XML_Retriever', level=log_level)

class UniProt_XML_Retriever(XML_Retriever):

    def __init__(self, retry_num=3, retry_wait=1, log_level='WARNING'):
        self.base_url = "https://www.uniprot.org/uniprot"
        self.retry_num = retry_num
        self.retry_wait = retry_wait
        self.logger = add_logger('UniRef_XML_Retriever', level=log_level)




