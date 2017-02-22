APP = imicrobe-blast-0.0.1
EMAIL = kyclark@email.arizona.edu

clean:
	find . \( -name \*.conf -o -name \*.out -o -name \*.log -o -name \*.param -o -name launcher-\* \) -exec rm {} \;

upload:
	files-upload -F $(WORK)/$(APP) kyclark/applications

update:
	apps-addupdate -F app.json

test: clean
	sbatch test.sh
