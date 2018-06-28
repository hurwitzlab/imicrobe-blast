APP = imicrobe-blast
SINGULARITY = singularity
SIZE = 1024
VERSION = 0.0.5
IMAGE_NAME = $(APP)-$(VERSION).img

img:
	rm -f $(IMAGE_NAME)
	sudo $(SINGULARITY) create --size $(SIZE) ~/$(IMAGE_NAME)
	sudo $(SINGULARITY) bootstrap ~/$(IMAGE_NAME) image.def

run:
	$(SINGULARITY) run $(IMAGE_NAME)

shell:
	sudo $(SINGULARITY) shell --writable -s /bin/bash $(IMAGE_NAME) 
	#sudo $(SINGULARITY) shell --writable -B $(PWD):/data -s /bin/bash $(IMAGE_NAME) 

scp:
	scp $(IMAGE_NAME) stampede:/work/03137/kyclark/$(APP)/stampede

iput:
	iput -f ~/$(IMAGE_NAME)

iget:
	cd /work/05066/imicrobe/singularity/; iget -f $(IMAGE_NAME); chmod ag+r ${IMAGE_NAME}
	irm $(IMAGE_NAME)