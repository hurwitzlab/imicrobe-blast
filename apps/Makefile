APP = imicrobe-blast
VERSION = 0.0.5
EMAIL = kyclark@email.arizona.edu

test:
	cd stampede2; sbatch test.sh

submit-test-job:
	jobs-submit -F stampede2/job.json

submit-test-job-to-public-app:
	jobs-submit -F stampede2/job-public-app.json

files-delete:
	files-delete $(CYVERSEUSERNAME)/applications/$(APP)-$(VERSION)

files-upload:
	files-upload -F imicrobe-blast/stampede2/run.sh $(CYVERSEUSERNAME)/applications/$(APP)-$(VERSION)
	files-upload -F imicrobe-blast/stampede2/template.sh $(CYVERSEUSERNAME)/applications/$(APP)-$(VERSION)

apps-addupdate:
	apps-addupdate -F imicrobe-blast/stampede2/app.json

deploy-app: files-delete files-upload apps-addupdate
