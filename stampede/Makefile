APP = muscope-blast-0.0.2
EMAIL = kyclark@email.arizona.edu

clean:
	find . \( -name \*.out -o -name \*.param \) -exec rm {} \;

upload: clean
	files-upload -F $(WORK)/$(APP) kyclark/applications

update:
	apps-addupdate -F app.json

test: clean
	sbatch test.sh

job:
	jobs-submit -F job.json

scriptsgz:
	(cd scripts && tar cvf ../scripts.tgz *)

jobs-template:
	jobs-template $(APP) > job.json

jobs-submit:
	jobs-submit -F job.json
