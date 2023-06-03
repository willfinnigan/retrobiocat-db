import pytest
from retrobiocat_web.app.app import create_app
from retrobiocat_web.app.retrobiocat.forms.retrosynthesis_options import Retrosynthesis_Config_Form


@pytest.fixture()
def app():
    app = create_app()
    app.config.update({
        "TESTING": True,
    })

    # other setup can go here

    yield app

    # clean up / reset resources here

@pytest.fixture()
def client(app):
    return app.test_client()

@pytest.fixture()
def runner(app):
    return app.test_cli_runner()

class Test_Retrosynthesis_Config_Form():

    def test_form_is_created_ok(self, client):
        ret_form = Retrosynthesis_Config_Form()
