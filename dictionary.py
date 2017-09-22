from Bio import SeqIO
import numpy as np

#Instanciacao das estruturas de dados do tipo dictionary
listOfSeq = {}
matrix = {}

p_value = {
		"ANDR" : "6.33741", "AP2A" : "5.40371", "AP2C" : "4.75071", "ASCL1" : "5.00301", "ATF3" : "4.59231", "ATF4" : "4.25481", "BACH2" : "5.36695", "BATF" : "6.12211", "BC11A" : "6.53131", "BCL6" : "6.24331", "BHE40" : "4.44851", "BMAL1" : "4.81285", "BRAC" : "5.93831", "CDX2" : "6.44232", "CEBPA" : "6.15857", "CEBPB" : "5.66791", "COE1" : "5.21891", "COT2" : "6.04551", "CREB1" : "4.96052", "CTCF" : "1.83531", "CTCFL" : "2.54771", "DUX4" : "4.24667", "E2F1" : "3.95651", "E2F3" : "5.60783", "E2F4" : "3.77691", "E2F6" : "4.97661", "EGR1" : "4.59371", "EGR2" : "4.65681", "ELF1" : "4.57601", "ELF3" : "5.55821", "ELF5" : "5.32491", "ELK4" : "4.08196", "ERG" : "5.52901", "ERR1" : "5.23021", "ESR1" : "4.76651", "ESR2" : "4.90741", "ETS1" : "5.12281", "ETV1" : "4.57559", "EVX2" : "6.16703", "FLI1" : "5.25931", "FOS" : "5.79736", "FOSB" : "5.84039", "FOSL1" : "4.61924", "FOSL2" : "4.66692", "FOXA1" : "6.01702", "FOXA2" : "5.81622", "FOXH1" : "5.95429", "FOXK1" : "6.5611", "FOXM1" : "6.47251", "FOXO1" : "6.44192", "FOXP1" : "6.51129", "GABPA" : "4.03691", "GATA1" : "5.73141", "GATA2" : "5.67971", "GATA3" : "6.65482", "GATA4" : "5.99156", "GATA6" : "5.73231", "GCR" : "5.45461", "GFI1B" : "6.34433", "GRHL2" : "5.53284", "HNF1B" : "5.56241", "HNF4A" : "5.36451", "HSF1" : "4.94931", "HTF4" : "5.76816", "HXB13" : "6.53543", "IRF1" : "5.29681", "IRF2" : "4.43121", "IRF4" : "5.19551", "ISL1" : "6.21911", "JUN" : "5.22518", "JUNB" : "5.65817", "JUND" : "5.26441", "KAISO" : "2.83851", "KLF1" : "3.56351", "KLF15" : "3.31911", "KLF4" : "4.0354", "KLF5" : "4.73701", "KLF6" : "4.62781", "LEF1" : "6.72381", "LHX2" : "6.62512", "LYL1" : "5.94061", "MAF" : "6.09041", "MAFG" : "2.37571", "MAFK" : "3.97231", "MAX" : "4.64031", "MAZ" : "4.64721", "MEF2A" : "6.47571", "MEF2B" : "6.59741", "MEF2C" : "6.42561", "MEF2D" : "5.99668", "MEIS1" : "6.27105", "MITF" : "4.00992", "MXI1" : "4.26221", "MYB" : "5.98733", "MYC" : "4.44066", "MYCN" : "4.20051", "MYOD1" : "4.47831", 
		"NANOG" : "6.62011", "NDF1" : "5.05856", "NF2L2" : "3.74991", "NFE2" : "4.69521", "NFIC" : "5.60251", "NFKB1" : "6.85701", "NFYA" : "5.39301", "NFYB" : "4.89731", "NFYC" : "5.90781", "NKX21" : "6.21355", "NR4A1" : "6.34763", "NRF1" : "0.44511", "OTX2" : "6.19477", "P53" : "3.44591", "P63" : "3.98641", "P73" : "5.15901", "PAX5" : "3.77741", "PBX1" : "5.72289", "PBX3" : "5.52922", "PDX1" : "6.37285", "PO2F2" : "6.29261", "PO3F2" : "6.32701", "PO5F1" : "6.18091", "PPARG" : "5.79441", "PRD14" : "5.47481", "PRDM1" : "5.90761", "PRGR" : "5.84971", "RARA" : "5.31441", "REST" : "0.74791", "RFX2" : "2.32281", "RFX5" : "4.95171", "RUNX1" : "5.67151", "RUNX2" : "5.91871", "RUNX3" : "5.54829", "RXRA" : "5.86461", "SIX1" : "5.83001", "SIX2" : "6.24631", "SMAD2" : "5.71925", "SNAI2" : "4.04502", "SOX2" : "6.70601", "SP1" : "-1.56869", "SP2" : "2.25691", "SP4" : "4.49461", "SPI1" : "4.93661", "SPIB" : "4.20911", "SRBP1" : "5.50945", "SRF" : "5.82391", "STA5A" : "5.80437", "STA5B" : "6.37151", "STAT1" : "5.26601", "STAT2" : "5.46721", "STAT3" : "6.02026", "STAT4" : "6.68682", "SUH" : "5.85721", "TAF1" : "5.02901", "TAL1" : "5.00481", "TBP" : "6.21194", "TBX21" : "5.78941", "TCF7" : "6.49011", "TEAD1" : "6.13019", "TEAD4" : "6.32661", "TF65" : "5.42531", "TF7L2" : "6.29291", "TFAP4" : "4.89738", "TFE2" : "4.94387", "TGIF1" : "5.8232", "TWST1" : "5.68941", "TYY1" : "3.52073", "USF1" : "2.73704", "USF2" : "2.05731", "VDR" : "5.19001", "ZBT17" : "5.44691", "ZBT7A" : "4.85817", "ZEB1" : "5.92264", "ZFP42" : "5.52841", "ZFX" : "4.84128", "ZN143" : "2.57551", "ZN263" : "4.12761", "ZN274" : "4.85761", "ZN281" : "3.37831", "ZN335" : "4.42511", "HSF2" : "5.80351", "FOXJ3" : "5.01711", "BACH1" : "4.46071", "CEBPE" : "6.23781", "ERR2" : "5.54155", "ARNT" : "3.78077", "ATF1" : "5.04362", "ATF2" : "5.40369", "ATOH1" : "5.19272", "CEBPG" : "5.18655", "E2F7" : "5.12342", "EHF" : "4.84231", "ELK1" : "4.37933", "EPAS1" : "4.18377", "ETV4" : "6.46311", "FEV" : "5.51391", "FOXO3" : "6.15962", "HNF4G" : "5.66891", "HXA9" : "6.31196", "HXB4" : "6.65045", "IRF3" : "6.24891", "KLF3" : "3.28141", 
		"MAFB" : "6.08586", "MAFF" : "5.88221", "MZF1" : "6.39972", "NDF2" : "5.53264", "NFAC1" : "6.82671", "NFKB2" : "5.45987", "NKX25" : "5.98306", "NKX61" : "6.72621", "NR1D1" : "6.01111", "NR1H3" : "5.76761", "NR2C2" : "5.12743", "NR5A2" : "5.53383", "PPARA" : "5.94401", "RARG" : "5.50731", "RFX1" : "2.53751", "SALL4" : "5.73848", "SMAD3" : "5.96431", "SMAD4" : "5.42291", "SOX10" : "6.42631", "SOX3" : "6.53663", "SOX4" : "6.48933", "SOX9" : "6.03181", "STAT6" : "6.34923", "STF1" : "5.81942", "TF7L1" : "6.66001", "TFE3" : "4.17543", "THA11" : "2.68361", "ZN322" : "3.96011", "DBP" : "5.98585", "MEIS2" : "4.40467", "EVI1" : "5.34161", "SRBP2" : "5.45601", "ZIC1" : "5.44086", "E2F5" : "4.42503", "NFAC3" : "6.43757", "SP3" : "3.73391", "E2F2" : "4.31021", "MBD2" : "2.75724", "ALX1" : "6.36661", "ATF6A" : "3.94011", "SRY" : "6.57201", "AP2B" : "4.9328", "AHR" : "3.99008", "FOXI1" : "4.76153", "ERR3" : "5.86695", "FOXA3" : "6.13292", "CRX" : "6.18351", "BHA15" : "5.08805", "NR6A1" : "5.35025", "HNF6" : "5.40454", "GLI3" : "5.60046", "MYOG" : "4.59741", "PTF1A" : "5.32611", "NR1H4" : "4.42431", "ZKSC1" : "4.00721", "ZIC3" : "4.82201", "ETV2" : "4.85041", "IRF8" : "5.03961", "OLIG2" : "6.08881", "BATF3" : "6.28661", "ETS2" : "6.57261", "REL" : "6.21151", "NFAC2" : "6.95952", "PKNX1" : "5.44108", "RFX3" : "2.61281", "RXRG" : "3.34671", "CEBPD" : "6.04011", "CLOCK" : "5.32721", "COT1" : "5.41421", "ELF2" : "4.77111", "ETV5" : "6.31411", "FEZF1" : "6.70421", "FOXP2" : "6.22857", "HIF1A" : "4.90321", "HNF1A" : "6.18801", "HXC9" : "6.50914", "ITF2" : "5.79733", "KLF12" : "3.12986", "KLF9" : "2.08171", "NKX31" : "6.37241", "OSR2" : "5.47301", "OZF" : "-8.22369", "PATZ1" : "4.32551", "PAX6" : "5.90291", "PO2F1" : "6.44341", "PRDM6" : "7.16031", "RORA" : "6.17471", "SMCA1" : "6.65377", "SMCA5" : "5.53681", "SOX17" : "6.38299", "TBX3" : "6.10511", "TFDP1" : "4.59911", "THAP1" : "3.15301", "VEZF1" : "5.22871", "Z324A" : "5.39071", "Z354A" : "6.16091", "ZBT14" : "2.16434", "ZBT18" : "5.68945", "ZBT48" : "5.85594", "ZBTB6" : "4.87951", 
		"ZFP28" : "4.39431", "ZFP82" : "5.87621", "ZIM3" : "4.41901", "ZN121" : "1.34901", "ZN134" : "4.45691", "ZN136" : "0.70801", "ZN140" : "2.89611", "ZN214" : "6.26261", "ZN250" : "-2.25579", "ZN257" : "5.31539", "ZN260" : "-0.86839", "ZN264" : "5.29771", "ZN317" : "4.96011", "ZN320" : "3.24141", "ZN329" : "-0.92689", "ZN331" : "0.89351", "ZN341" : "4.92701", "ZN350" : "6.77721", "ZN382" : "3.57251", "ZN394" : "6.90731", "ZN418" : "6.02501", "ZN436" : "-4.20249", "ZN449" : "5.25821", "ZN467" : "3.41481", "ZN490" : "-5.50589", "ZN502" : "0.15941", "ZN528" : "3.45451", "ZN547" : "4.47421", "ZN549" : "4.40461", "ZN554" : "5.16471", "ZN563" : "5.24111", "ZN582" : "3.13841", "ZN586" : "-4.80379", "ZN667" : "5.66511", "ZN680" : "1.56241", "ZN708" : "3.98121", "ZN768" : "2.36371", "ZN770" : "3.07541", "ZN816" : "4.10191", "ZNF18" : "4.99364", "ZNF41" : "4.35271", "ZNF76" : "3.57581", "ZNF8" : "-1.55239", "ZNF85" : "5.47241", "ZSC22" : "3.90511", "ZSC31" : "4.62081", "NFIA" : "5.75961", "NR2C1" : "5.63212", "TFEB" : "1.52516", "PIT1" : "6.63531", "AIRE" : "6.37171", "NR1I2" : "6.04901", "NR1I3" : "5.80501", "MTF1" : "3.76681", "HXA13" : "6.12954", "IKZF1" : "6.3404", "ARI5B" : "5.72941", "HXA10" : "6.17423", "CREM" : "5.02699", "HIC1" : "5.36472", "NF2L1" : "6.11632", "NFAC4" : "6.72218", "RELB" : "5.93401", "NR4A2" : "6.15657", "CDX1" : "5.93803", "CUX1" : "5.52101", "FOXO4" : "6.50759", "NOBOX" : "6.46625", "ZN384" : "6.771", "HEN1" : "2.95401", "DLX3" : "6.13416", "MYF6" : "5.47889", "OVOL1" : "5.45742", "FOXJ2" : "6.54886", "PEBB" : "5.82942", "PRRX2" : "6.88174", "IRF9" : "6.26288", "SOX5" : "6.32764", "RXRB" : "5.73829", "HXB8" : "6.48744", "FOXC1" : "5.93541", "GFI1" : "6.43254", "HINFP" : "3.36791", "MECP2" : "3.25632", "NR2E3" : "6.55341", "PBX2" : "6.44341", "HXA1" : "6.09561", "KLF8" : "5.29075", "IRF7" : "6.62174", "SNAI1" : "5.33066", "NKX28" : "5.88086", "FOXQ1" : "6.38518", "INSM1" : "4.81366", "HXB7" : "6.5474", "THB" : "5.27371", "WT1" : "2.91481", "THA" : "5.66091", "NKX32" : "5.95692", "HLF" : "5.35061", "LHX3" : "6.23335", "PO3F1" : "6.33031", "RORG" : "5.93884"
	}	

#Leitura do arquivo fasta e separacao das sequencias validas	
def readFastaFile(fileIn):
	for seq_record in SeqIO.parse(fileIn, "fasta"):
		if (len(seq_record) > 0):
			key = str(seq_record.id)
			sequence = str(seq_record.seq)
			listOfSeq.update({key : sequence.upper()})

#Adiciona aa arvore matrix o valor "value" na chave "keyName"
def addToMatrix(keyName, value):
	value = value.replace('\t', ' ')
	value = value.replace('\n', '; ')
	value = value[:-2] #Utilado para remover o '\n' da ultima linha. 2 porque '\n' eh substituido por '; ', 2 characteres na ultima linha
	value = np.matrix(value)
	value = np.round(value, decimals=3)
	matrix[keyName] = value


#Le um arquivo no formato de matrizes do site hocomoco
def readMatrix(fileIn):
	with open(fileIn, "r") as matrixFile:
		i = 0
		matriz = ''
		name = ''
		for line in matrixFile:
			if line.startswith('>'):
				if i>0:
					addToMatrix(name, matriz)
				name, sep, lixo = line.partition('_')
				name = name.lstrip('>')
				matriz=''
				i=i+1
			else:
				matriz += line
		addToMatrix(name, matriz) #adding last one			

def row(matrix, index):
   	return v_Matx[index]

def getLetterValue(letter, array):
	if letter == 'A':
		return array[0]
	elif letter == 'C':
		return array[1]
	elif letter == 'G':
		return array[2]
	elif letter == 'T':
		return array[3]	
	else:
		print('ERROR: I dont know how to compare letter ' + letter)
		#return 0 #How to treat it?


def compareMotif(motif, matrixMotif):
	motifValue = 0
	for i in range(0, len(motif)):
		motifValue+=getLetterValue(motif[i], row(matrixMotif, i))
	return motifValue

def findAndCompareMotif(sequence, matrixMotif, k_Matx):
	motif = ''
	toOut=''
	#Finding Motif
	for i in range(0, len(sequence)-len(matrixMotif)+1):	
		motif = sequence[i : i+len(matrixMotif)]
		aux=compareMotif(motif, matrixMotif)
			
		if aux >= p_value[k_Matx]: 
			toOut+=str(aux)+'; '
			print(aux)
			print('>=')
			print(p_value[k_Matx])
		motif='' #Cleaning Motif
	return toOut


#Funcao principal
if __name__ == '__main__':
	readFastaFile("./fasta/pgene.txt")
	readMatrix("./matrixMotif/pwms_HUMAN_mono.txt")
	i=1
	out = open('out.txt', 'w')
	for k_Seq,v_Seq in listOfSeq.items():
		print(str(i)+' '+str(k_Seq))
		i+=1
		for k_Matx,v_Matx in matrix.items():
			motif = findAndCompareMotif(v_Seq, v_Matx, k_Matx)
			if motif != '':
				out.write(' > '+k_Matx+' > '+motif[:-2])
				out.write('\n')
			
	out.close()			
	



	#printListOfSeq()