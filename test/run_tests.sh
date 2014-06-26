#!/bin/sh

CLUSTER=../cluster
TAGS=("ACTG" "CAGT" "TGAC" "GTCA")

FAILURES=0

function check_tag_parse {
	TAGFILE=$1

	echo "Checking output with tagfile ${TAGFILE}"

	${CLUSTER} tag_parse ${TAGFILE} AAAA GGGG 2> cluster_output.txt
	grep " using 4 tags for this FASTQ." cluster_output.txt > /dev/null
	if [[ $? -ne 0 ]];
	then
		let "FAILURES=FAILURES+1"
		echo "FAILURE: tags didn't load properly."
	fi

	function expect_file_length {
		FILENAME=$1
		EXPECTED_LINES=$2
		LINES=$(wc -l < $1)

		if [[ ${LINES} -ne ${EXPECTED_LINES} ]] ;
		then 
			let "FAILURES=FAILURES+1"
			echo "FAILURE: found wrong number of lines in ${FILENAME} file"
		fi
	}

	for tag in "${TAGS[@]}"
	do
		expect_file_length "tag_parse.${tag}.txt" 1
		expect_file_length "tag_parse.${tag}_clusters.csv" 1
	done

	rm tag_parse*.txt *.csv
	rm cluster_output.txt
}

#first, test to ensure that tag files are being parsed
check_tag_parse tag_parse.tags
check_tag_parse tag_parse2.tags
check_tag_parse tag_parse3.tags
check_tag_parse tag_parse4.tags

if [[ ${FAILURES} -ne 0 ]];
then
	echo "Tests failed (${FAILURES} failures)"
else
	echo "Tests passed"
fi
