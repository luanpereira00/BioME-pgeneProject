import time
from thread import *
from concurrent.futures import *

a_lock = allocate_lock()
executor = ThreadPoolExecutor(max_workers=6)

def funcaoTeste():
	with a_lock:		
		print("a_lock is locked while this executes")
		print("Sou uma funcao de thread")
		time.sleep(2)



#Funcao principal
if __name__ == '__main__':
	print("Estou na funcao principal")
	a = executor.submit(funcaoTeste)
	a.result(timeout=1)
	b = executor.submit(funcaoTeste)
	c = executor.submit(funcaoTeste)
	d = executor.submit(funcaoTeste)
	e = executor.submit(funcaoTeste)
	f = executor.submit(funcaoTeste)
	g = executor.submit(funcaoTeste)
	h = executor.submit(funcaoTeste)
	i = executor.submit(funcaoTeste)
