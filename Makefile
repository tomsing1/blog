all: render push

render:
	quarto render

push:
	git push
