
from os import listdir
from os import remove
from sys import argv
from os import path
import subprocess
import pyperclip

from colorama import Fore, Style, init

init()
print(Style.BRIGHT, end='')
if len(argv) < 2 or len(argv) > 3:
	print(Fore.RED + "ERROR" + Fore.RESET)
	print("Usage: run.py <executable> [input #]")
	exit(0)

for filename in listdir():
	if filename.endswith('.off') or filename == "MST.txt":
		remove(filename)

inputFiles = []
with open('input.txt') as file:
	for line in file:
		inputFiles.append([line.split()[0], line.split()[3]])

choice = 0

if len(argv) == 2:
	print("Choose input File to run", Fore.RED, argv[1], Fore.RESET, "with")
	for i in range(len(inputFiles)):
		if inputFiles[i][1] == "total":
			print((Fore.GREEN + '0.' + Fore.RESET + ' EXIT'))
		else:
			print((Fore.GREEN + '{:2d}.' + Fore.RESET + ' {:7s} points {:s}').format(\
                            i + 1, inputFiles[i][0], inputFiles[i][1]))
	choice = int(input())
else:
	choice = int(argv[2])

if choice == 0:
	exit(0)

executable = argv[1].replace('/', '\\')
outputFile = 'output_' + path.basename(argv[1]).replace(".exe", "") + '.off'
logFile = 'log_' + path.basename(argv[1]).replace(".exe", "") + '.txt'
inputFile = inputFiles[choice - 1][1].replace('/', '\\')

print(executable, inputFile, outputFile, ">", logFile)
pyperclip.copy(executable + ' ' + inputFile + ' ' + outputFile + ' ' + ">" +
               ' ' + logFile)

print(Style.RESET_ALL)
