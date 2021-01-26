from requests import Session
from os.path import isdir, expanduser
from getpass import getpass
import xml.etree.ElementTree as ET


class CERNSession(Session):
    LOGIN_URL = "https://login.cern.ch"
    CERT_DIRS = [
        "/etc/grid-security/certificates",
        "/cvmfs/grid.cern.ch/etc/grid-security/certificates",
    ]

    def __init__(
        self,
        key_file="~/.globus/userkey.pem",
        cert_file="~/.globus/usercert.pem",
        ca_cert_dir=True,
    ):
        super().__init__()

        if ca_cert_dir is True:
            for dir in filter(isdir, self.CERT_DIRS):
                ca_cert_dir = dir
                break
            else:
                ca_cert_dir = None

        c = self.get_adapter(self.LOGIN_URL).get_connection(self.LOGIN_URL)
        c.cert_file = expanduser(cert_file)
        c.key_file = expanduser(key_file)
        if ca_cert_dir:
            c.ca_cert_dir = expanduser(ca_cert_dir)
        c.key_password = lambda: getpass("Password (%s): " % key_file)
        assert set(c.pool.queue) == {None}

        self.headers["User-Agent"] = "curl-sso-certificate/0.6"

    def _get_form(self, resp):
        if resp.status_code == 200 and resp.url.startswith(self.LOGIN_URL):
            return ET.fromstring(resp.content).find("body/form")
        else:
            return None

    def get_redirect_target(self, resp):
        url = super().get_redirect_target(resp)
        if url is None:
            form = self._get_form(resp)
            if form:
                url = form.get("action")
        return url

    def rebuild_auth(self, preq, resp):
        super().rebuild_auth(preq, resp)
        form = self._get_form(resp)
        if form:
            preq.prepare_method("POST")
            preq.prepare_body({
                el.get("name"): el.get("value")
                for el in form.findall("input")
            }, None)
