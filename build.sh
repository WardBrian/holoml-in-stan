#!/bin/sh
jupyter nbconvert "./src/HoloML in Stan.ipynb" --to html --output-dir="./docs" --template classic --TagRemovePreprocessor.remove_input_tags=hide-code -CSSHTMLHeaderPreprocessor.style=tango --execute
