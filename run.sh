#! /bin/bash

red=$'\e[1;31m'
grn=$'\e[1;32m'
yel=$'\e[1;33m'
blu=$'\e[1;34m'
mag=$'\e[1;35m'
cyn=$'\e[1;36m'
end=$'\e[0m'

rm tempMod.txt*

# if [ 1 -eq $# ]

if [ "$OSTYPE" != "msys" ]
then
	make $1 
	if [ $? -ne 0 ]
	then
    	echo Make Failed
    	exit 1
	fi
fi


if [ 1 -eq $# ]
then
let i=0
echo "Chose input File to run ${red}$1${end} with" 
while IFS=' ' read -u 10 -a line; do
	let i++
	if [ ${line[3]} == "total" ] 
	then
		printf "${grn}%2d.${end} EXIT" 0
	else
		printf "${grn}%2d.${end} %10d points %10s" $i ${line} ${line[3]}
	fi
	echo
done 10<input.txt
	read choice
	choice=$(($choice+0))
else
	choice=$(($2+0))
fi
if [ 0 -eq $choice ] 
then
	exit 0
fi
let i=0
fileName=""
while IFS=' ' read -u 10 -a line; do
	let i++
	if [ $i -eq $choice ] 
	then
		fileName=${line[3]}
		break
	fi
done 10<input.txt
printf "${red}Running Command${end}"
echo
printf "${mag}./$1 $fileName output_"$(basename $1 .exe)".txt > log_"$(basename $1 .exe)".txt${end}"
echo
if [[ "$OSTYPE" == "msys" ]]
then
	./$1 $fileName output_"$(basename $1 .exe)".txt > log_"$(basename $1 .exe)".txt
	tail -n1 log_"$(basename $1 .exe)".txt
else
	./$1 $fileName output_$1.txt > log_$1.txt
	tail -n1 log_$1.txt
fi
