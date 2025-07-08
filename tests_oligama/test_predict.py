from ollikit.oligama.predictors.affinity_predictors import Nupack_Affinity_Predictor

# Пример двух последовательностей (строки)
seq1 = "ACTGCTAGAGATTTTCCACAT"
seq2 = "TGTGGAAAATCTCTAGCAGTT"

# Инициализация предиктора
predictor = Nupack_Affinity_Predictor(celsius=25, material='dna')

# Предсказание аффинности для двух строк
result = predictor.predict(seq1, seq2, units='fraction')
print("Результат аффинности для двух строк:", result)

# Пример с массивами строк
seqs1 = ["ACTGCTAGAGATTTTCCACAT", "AATCGCTAGCTAGCTAGCTAG"]
seqs2 = ["TGTGGAAAATCTCTAGCAGTT", "ACTGCTAGAGATTTTCCACAT"]

result_arr = predictor.predict(seqs1, seqs2, units='fraction')
print("Результаты аффинности для массивов строк:", result_arr) 