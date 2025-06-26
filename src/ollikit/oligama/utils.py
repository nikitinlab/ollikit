import itertools
import numpy as np
import random
import keras
import tensorflow as tf
import tensorflow.keras.saving
import tensorflow.keras.backend as K

import pandas as pd



def df_to_excel(df_arr, sheet_names, path, scientific_format_flag = False):
	"""
	Сохраняет список DataFrame в файл Excel.

	Args:
		df_arr (list): Список DataFrame.
		sheet_names (list): Список имен листов.
		path (str): Путь к файлу Excel.
		scientific_format_flag (bool): Флаг для использования научного формата чисел.
	"""    
	with pd.ExcelWriter(path, 
					  	engine='xlsxwriter',
						mode='w') as writer:
		if scientific_format_flag:
			format1 = writer.book.add_format({'num_format': '0.00E+00'})
		else:
			format1 = writer.book.add_format({'num_format': '0.00'})
		for i, df in enumerate(df_arr):
			df.to_excel(writer, sheet_name=sheet_names[i], index = False)
			for column in df.columns:
				column_length = max(df[column].astype(str).map(len).max(), len(column)) + 5
				col_idx = df.columns.get_loc(column)
				writer.sheets[sheet_names[i]].set_column(col_idx, col_idx, column_length, format1)



@keras.saving.register_keras_serializable() # Для сериализации/десериализации модели
def custom_loss(y_true, y_pred):
	"""
	Пользовательская функция потерь.

	Args:
		y_true: Истинные значения.
		y_pred: Предсказанные значения.

	Returns:
		tf.Tensor: Значение потерь.
	"""    
	y_true = tf.where(tf.math.is_nan(y_true), y_pred, y_true)

	mse_loss = K.mean(tf.square(y_true - y_pred), axis=-1) 
	constraint_loss =  K.mean(tf.maximum(0.0, (y_pred[:, 1] - y_pred[:, 0])))
       
	mse_loss = K.cast(mse_loss, "float32")
	constraint_loss = K.cast(constraint_loss, "float32")
	return mse_loss + constraint_loss

def init_kmer_vocab(kmer_len, material='dna'):
	"""
    Создает словарь k-мерных пар.

    Args:
        kmer_len (int): Длина k-мера.
        material (str): Тип материала ('dna' или 'rna').

    Returns:
        list: Список уникальных k-мерных пар.
	"""    
	letters = ['A', 'T', 'G', 'C'] if material.lower() == 'dna' else ['A', 'U', 'G', 'C']
	kmer_pairs = ["".join(x) for x in list((itertools.product(*[letters]*kmer_len*2)))]
	unique_kmer_pairs = []
	for kmer_pair in kmer_pairs:
		kmer1 = kmer_pair[:len(kmer_pair)//2]
		kmer2 = kmer_pair[len(kmer_pair)//2:]
		synonyms = [kmer1+kmer2, kmer2+kmer1, (kmer1+kmer2)[::-1], (kmer2+kmer1)[::-1]]
		if np.array([syn in unique_kmer_pairs for syn in synonyms]).any():
			continue
		unique_kmer_pairs.append(kmer1 + kmer2)
	return unique_kmer_pairs

def encode_kmers(seq1, seq2, unique_kmer_pairs):
	"""
    Кодирует последовательности в k-мерные пары.

    Args:
        seq1 (str): Первая последовательность.
        seq2 (str): Вторая последовательность.
        unique_kmer_pairs (list): Список уникальных k-мерных пар.

    Returns:
        str: Закодированная строка.
	"""
	kmer_len = len(unique_kmer_pairs[0])//2
	seq2 = seq2[::-1]
	encoded = []
	for i in range(len(seq1) - kmer_len + 1):
		kmer_pair = seq1[i:i+kmer_len] + seq2[i:i+kmer_len]

		kmer1 = kmer_pair[:len(kmer_pair)//2]
		kmer2 = kmer_pair[len(kmer_pair)//2:]
		synonyms = [kmer1+kmer2, kmer2+kmer1, (kmer1+kmer2)[::-1], (kmer2+kmer1)[::-1]]
		for syn in synonyms:
			if syn in unique_kmer_pairs:
				encoded.append(syn)
				break
	return " ".join(encoded)

def random_seq(length, material='dna'):
	"""
    Генерирует случайную последовательность ДНК.

    Args:
        length (int): Длина последовательности.
 		material (str): Тип материала ('dna' или 'rna').
   
    Returns:
        str: Случайная последовательность.
	"""        
	letters = ['A', 'T', 'G', 'C'] if material.lower() == 'dna' else ['A', 'U', 'G', 'C']
	return "".join(random.choices(letters, k=length))

def compl(seq, material='dna'):
	"""
    Возвращает комплементарную последовательность ДНК или РНК.

    Args:
        seq (str): Последовательность ДНК или РНК.
        material (str): Тип материала ('dna' или 'rna').

    Returns:
        str: Комплементарная последовательность.
	"""    

	compl_dict = {"A": "T", "T": "A", "G": "C", "C": "G"} if material.lower() == 'dna' else {"A": "U", "U": "A", "G": "C", "C": "G"}
	return "".join(compl_dict[let] for let in seq)[::-1]



def mutate_letter(seq, ind, material='dna'):
	"""
    Заменяет букву в последовательности на другую случайную букву.

    Args:
        seq (str): Последовательность.
        ind (int): Индекс буквы для замены.
        material (str): Тип материала ('dna' или 'rna').

    Returns:
        str: Измененная последовательность.
	"""    
	letters = ['A', 'T', 'G', 'C'] if material.lower() == 'dna' else ['A', 'U', 'G', 'C']
	random.seed()
 
	# Проверяем, что символ принадлежит допустимому набору
	#if seq[ind] not in letters:
	#	raise ValueError(f"Invalid symbol '{seq[ind]}' in sequence at index {ind} for material '{material}'. Expected one of {letters}.")
 
	letters.remove(seq[ind])
	mod_let = letters[random.randint(0, 2)]
	mod_seq = list(seq)
	mod_seq[ind] = mod_let
	
	return "".join(mod_seq)

def mutate_x_letters(seq, x, material='dna'):
    """
    Вносит x случайных мутаций в последовательность.

    Args:
        seq (str): Последовательность.
        x (int): Количество мутаций.
        material (str): Тип материала ('dna' или 'rna').

    Returns:
        str: Измененная последовательность.
    """    
    random.seed()
    # Ограничиваем x длиной последовательности
    x = min(x, len(seq))
    mut_list = random.sample(range(0, len(seq)), x)
    mod_seq = seq
    for ind in mut_list:
        mod_seq = mutate_letter(mod_seq, ind, material)
    return mod_seq


def crossover(donor_seq, acceptor_seq, max_batch_size = None, num_cycles = 1):
	"""
    Выполняет скрещивание двух последовательностей.

    Args:
        donor_seq (str): Донорская последовательность.
        acceptor_seq (str): Акцепторная последовательность.
        max_batch_size (int): Максимальный размер батча для скрещивания.
        num_cycles (int): Количество циклов скрещивания.

    Returns:
        str: Результирующая последовательность.
	"""    
	if max_batch_size is None:
		max_batch_size = len(donor_seq)//2
	new_seq = acceptor_seq
	for i in range(num_cycles):
		batch_size = random.randint(1, max_batch_size)
		position = random.randint(0, len(donor_seq) - batch_size)
		new_seq = new_seq[:position] + donor_seq[position:position+batch_size] + new_seq[position+batch_size:]
	return new_seq

def multiple_crossover(seq_list, max_mutations=10, crossover_cycles=2, max_crossover_batch=8, material='dna'):
    """
    Выполняет множественное скрещивание и мутацию.

    Args:
        seq_list (list): Список последовательностей.
        max_mutations (int): Максимальное количество мутаций.
        crossover_cycles (int): Количество циклов скрещивания.
        max_crossover_batch (int): Максимальный размер батча для скрещивания.

    Returns:
        str: Результирующая последовательность.
    """

    new_seq = compl(seq_list[0])  # Создаём комплементарную последовательность
    for seq in seq_list[1:]:
        new_seq = crossover(new_seq, compl(seq), max_crossover_batch, crossover_cycles)
    
    # Ограничиваем количество мутаций длиной новой последовательности
    max_mutations = min(len(new_seq), max_mutations)
    return mutate_x_letters(new_seq, random.randint(1, max(1, max_mutations)), material)



def GC_content(seq):
	"""
    Вычисляет GC-состав последовательности.

    Args:
        seq (str): Последовательность ДНК или РНК.

    Returns:
        float: GC-состав.
	# """    
	# n = len(seq)
	# GC_count = 0
	# for let in seq:
	# 	if let in ['G', 'C']:
	# 		GC_count += 1
	# return GC_count/len(seq)
 
	return sum(1 for let in seq if let in "GC") / len(seq)