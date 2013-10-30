
def test_best():
	'''
	--- Doctest ---

	>>> import subprocess
	>>> p = subprocess.Popen(["python", "../scripts/main.py", "first_test/contigs.fa", "first_test/links.txt", "first_test/subseqs.txt"], stdout=subprocess.PIPE)
	>>> out, err = p.communicate()
	>>> out 
	'>scaffold1\\nAAAnCCCnGGGnTTTn\\n'
	>>> p = subprocess.Popen(["python", "../scripts/main.py", "second_test/contigs.fa", "second_test/links.txt", "second_test/subseqs.txt"], stdout=subprocess.PIPE)
	>>> out, err = p.communicate()
	>>> out 
	'>scaffold1\\nTTCnCAAn\\n'
	'''
	
	pass

if __name__ == "__main__":
    import doctest
    doctest.testmod()