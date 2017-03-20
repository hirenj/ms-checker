#!/bin/bash

dir_resolve()
{
cd "$1" 2>/dev/null || return $?  # cd to desired directory; if fail, quell any error messages but return exit status
echo "`pwd -P`" # output full, link-resolved path
}

optspec=":ion-:"
while getopts "$optspec" optchar; do
    case "${optchar}" in
        -)
            case "${OPTARG}" in
                input)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    input_dir="`dir_resolve \"$val\"`"
                    ;;
                input=*)
                    val=${OPTARG#*=}
                    opt=${OPTARG%=$val}
                    input_dir="`dir_resolve \"$val\"`"
                    ;;
                nonetwork)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    nonetwork="--offline"
                    ;;
                nonetwork=*)
                    val=${OPTARG#*=}
                    opt=${OPTARG%=$val}
                    nonetwork="--offline"
                    ;;
                output)
                    val="${!OPTIND}"; OPTIND=$(( $OPTIND + 1 ))
                    output_dir="`dir_resolve \"$val\"`"
                    ;;
                output=*)
                    val=${OPTARG#*=}
                    opt=${OPTARG%=$val}
                    output_dir="`dir_resolve \"$val\"`"
                    ;;
                *)
                    if [ "$OPTERR" = 1 ] && [ "${optspec:0:1}" != ":" ]; then
                        echo "Unknown option --${OPTARG}" >&2
                    fi
                    ;;
            esac;;
        h)
            echo "usage: $0 [--input[=]<value>] [--output[=]<value>]" >&2
            exit 2
            ;;
        n)
            nonetwork="--offline"
            ;;
        *)
            if [ "$OPTERR" != 1 ] || [ "${optspec:0:1}" = ":" ]; then
                echo "Non-option argument: '-${OPTARG}'" >&2
            fi
            ;;
    esac
done

if [[ -z "$input_dir" || -z "$output_dir" ]]; then
	echo "Could not find paths Input '$input_dir', Output '$output_dir'"
	exit 1
fi

if [ ! -z `which grunt` ]; then
    grunt version
fi

for recipe_dir in `find "${input_dir}" -type f -name '*.recipe.json' -exec dirname {} \; | grep -v '^\.' | uniq`; do
    runrecipe --input "$recipe_dir" --output "$output_dir" --nomangle
    if [ $? -gt 0 ]; then
        exit "$?"
    fi
done

find "${input_dir}" -name "manifest*.xlsx" -exec node js/main.js --manifest {} --outputdir "$output_dir" "$nonetwork" --prefix-basename \;


