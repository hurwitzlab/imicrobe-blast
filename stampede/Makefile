APP = imicrobe-blast-0.0.2
EMAIL = kyclark@email.arizona.edu

clean:
	find . \( -name \*.out -o -name \*.param \) -exec rm {} \;

up: 
	copy-from-manifest.pl6 $(WORK)/$(APP)

app:
	apps-addupdate -F app.json

test: clean
	sbatch test.sh

job:
	jobs-submit -F job.json

template:
	jobs-template $(APP) > job.json
