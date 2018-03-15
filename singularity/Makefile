APP = imicrobe-blast
SINGULARITY=singularity
SIZE = 1024
VERSION = 0.0.5
IMAGE_NAME = $(APP)-$(VERSION).img

img:
	rm -f $(IMAGE_NAME)
	sudo $(SINGULARITY) create --size $(SIZE) $(IMAGE_NAME)
	sudo $(SINGULARITY) bootstrap $(IMAGE_NAME) image.def

run:
	$(SINGULARITY) run $(IMAGE_NAME)

shell:
	sudo $(SINGULARITY) shell --writable -s /bin/bash $(IMAGE_NAME) 
	#sudo $(SINGULARITY) shell --writable -B $(PWD):/data -s /bin/bash $(IMAGE_NAME) 

scp:
	scp $(IMAGE_NAME) stampede:/work/03137/kyclark/$(APP)/stampede

iput:
	iput $(IMAGE_NAME)

iget:
	iget $(IMAGE_NAME)