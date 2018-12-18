#!/bin/sh


export opts='-shell-escape -pdf'
model='NEMO'

clean() {
    ## Delete latex build files
    find latex -regextype posix-extended                                              \
         -regex ".*\.(aux|bbl|blg|dvi|fdb|fls|idx|ilg|ind|log|maf|mtc|out|pdf|toc).*" \
         -exec rm {} \;

    ## Remove 'minted' directories
    find latex -type d -name '_minted*' -exec rm -r {} \;

    ## HTML exports
    find latex -type d -name 'html*'    -exec rm -r {} \;
}

build() {
    cd latex/$1/main
    latexmk $opts $1'_manual' > /dev/null
    mv            $1'_manual'.pdf ../../..
    cd -
}

clean
build $model

exit 0
