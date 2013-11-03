
def test_besst():
	'''
	--- Doctest ---

	>>> import subprocess
	>>> p = subprocess.Popen(["python", "../scripts/main.py", "first_test/contigs.fa", "first_test/links.txt", "first_test/subseqs.txt","0"], stdout=subprocess.PIPE)
	>>> out, err = p.communicate()
	>>> out 
	'>scaffold1\\nAAAnCCCnGGGnTTTn\\n'
	>>> p = subprocess.Popen(["python", "../scripts/main.py", "second_test/contigs.fa", "second_test/links.txt", "second_test/subseqs.txt", "0"], stdout=subprocess.PIPE)
	>>> out, err = p.communicate()
	>>> out 
	'>scaffold1\\nGAAnCAAn\\n'
	>>> p = subprocess.Popen(["python", "../scripts/main.py", "third_test/contigs.fa", "third_test/links.txt", "third_test/subseqs.txt", "4"], stdout=subprocess.PIPE)
	>>> out, err = p.communicate()
	>>> out 
	'>scaffold1\\nGAGAGANNNNCTCTCTnCTCTTCCAAn\\n'
	>>> p = subprocess.Popen(["python", "../scripts/main.py", "fourth_test/contigs.fa", "fourth_test/links.txt", "fourth_test/subseqs.txt", "0"], stdout=subprocess.PIPE)
	>>> out, err = p.communicate()
	>>> strings = out.split() 
	>>> '>scaffold1' in strings
	True
	>>> '>scaffold2' in strings
	True
	>>> 'CCNNCCNCTn' in strings
	True
	>>> 'AANNGTn' in strings
	True
	>>> p = subprocess.Popen(["python", "../scripts/main.py", "fifth_test/contigs.fa", "fifth_test/links.txt", "fifth_test/subseqs.txt", "3"], stdout=subprocess.PIPE)
	>>> out, err = p.communicate()
	>>> out 
	'>scaffold1\\nCTCTCTn\\n>scaffold2\\nGAGAGANNNNNCTCTTCCAAn\\n'
	>>> p = subprocess.Popen(["python", "../scripts/main.py", "fifth_test/contigs.fa", "fifth_test/links.txt", "fifth_test/subseqs.txt", "4"], stdout=subprocess.PIPE)
	>>> out, err = p.communicate()
	>>> out 
	'>scaffold1\\nGAGAGAnCTCTCTNNCTCTTCCAAn\\n'
	'''

	pass

if __name__ == "__main__":
    import doctest
    doctest.testmod()