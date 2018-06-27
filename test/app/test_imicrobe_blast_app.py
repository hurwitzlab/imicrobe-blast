import json

import pytest

from agave_app_test import submit_job, wait_for_job_to_finish


def test_imicrobe_blast_private(cyverse_oauth_session):
    job_json = json.loads("""\
        {
          "name": "imicrobe-blast-test",
          "appId": "imicrobe-blast-0.0.5",
          "defaultMaxRunTime": "00:30:00",
          "archive": true,
          "inputs": {
            "QUERY": "agave://data.iplantcollaborative.org/jklynch/data/muscope/last/test.fa"
          },
          "parameter": {
          }
        }
    """)

    submit_job_response = submit_job(
        cyverse_oauth_session=cyverse_oauth_session,
        job_json=job_json)
    print(submit_job_response.json())

    job_status = wait_for_job_to_finish(
        job_id=submit_job_response.json()['result']['id'],
        cyverse_oauth_session=cyverse_oauth_session)

    assert job_status == 'FINISHED'
