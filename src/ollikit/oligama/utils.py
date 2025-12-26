import itertools
import numpy as np
import random
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

    new_seq = compl(seq_list[0], material=material)  # Создаём комплементарную последовательность
    for seq in seq_list[1:]:
        new_seq = crossover(new_seq, compl(seq, material=material), max_crossover_batch, crossover_cycles)
    
    # Ограничиваем количество мутаций длиной новой последовательности
    max_mutations = min(len(new_seq), max_mutations)
    return mutate_x_letters(new_seq, random.randint(1, max(1, max_mutations)), material)


def _get_mutable_positions(template):
    """
    Возвращает список индексов позиций, которые можно изменять (x в шаблоне).
    
    Args:
        template: Шаблон последовательности
        
    Returns:
        list: Список индексов переменных позиций
    """
    return [i for i, char in enumerate(template) if char == 'x']


def _get_fixed_positions(template):
    """
    Возвращает словарь {индекс: символ} для фиксированных позиций.
    
    Args:
        template: Шаблон последовательности
        
    Returns:
        dict: Словарь фиксированных позиций
    """
    return {i: char for i, char in enumerate(template) if char != 'x' and char in 'ATGCU'}


def mutate_x_letters_with_template(seq, x, template="", material='dna'):
    """
    Вносит x случайных мутаций в последовательность с учетом шаблона.
    Не изменяет фиксированные позиции шаблона.
    
    Args:
        seq: Последовательность
        x: Количество мутаций
        template: Шаблон (если пустой, работает как обычная мутация)
        material: Тип материала
        
    Returns:
        str: Измененная последовательность
    """
    if not template:
        return mutate_x_letters(seq, x, material)
    
    # Получаем список позиций, которые можно мутировать
    mutable_positions = _get_mutable_positions(template)
    
    if not mutable_positions:
        # Все позиции фиксированы, мутации невозможны
        return seq
    
    # Ограничиваем x количеством доступных позиций
    x = min(x, len(mutable_positions))
    
    if x == 0:
        return seq
    
    # Выбираем случайные позиции только из мутабельных
    mut_list = random.sample(mutable_positions, x)
    
    mod_seq = list(seq.upper())
    letters = ['A', 'T', 'G', 'C'] if material.lower() == 'dna' else ['A', 'U', 'G', 'C']
    
    for ind in mut_list:
        if ind < len(mod_seq):
            # Убираем текущий символ из списка возможных
            available_letters = [l for l in letters if l != mod_seq[ind]]
            if available_letters:
                mod_seq[ind] = random.choice(available_letters)
    
    return "".join(mod_seq)


def crossover_with_template(donor_seq, acceptor_seq, template="", max_batch_size=None, num_cycles=1):
    """
    Выполняет скрещивание с учетом фиксированных позиций шаблона.
    
    Args:
        donor_seq: Донорская последовательность
        acceptor_seq: Акцепторская последовательность
        template: Шаблон
        max_batch_size: Максимальный размер батча
        num_cycles: Количество циклов
        
    Returns:
        str: Результирующая последовательность
    """
    if not template:
        return crossover(donor_seq, acceptor_seq, max_batch_size, num_cycles)
    
    if max_batch_size is None:
        max_batch_size = len(donor_seq) // 2
    
    new_seq = list(acceptor_seq.upper())
    fixed_positions = _get_fixed_positions(template)
    mutable_positions = _get_mutable_positions(template)
    
    # Выполняем скрещивание только в мутабельных позициях
    for i in range(num_cycles):
        if not mutable_positions:
            break  # Нет позиций для скрещивания
        
        batch_size = random.randint(1, min(max_batch_size, len(mutable_positions)))
        
        # Выбираем позицию только из мутабельных
        available_start_positions = [pos for pos in mutable_positions if pos + batch_size <= len(new_seq)]
        
        if not available_start_positions:
            break
        
        position = random.choice(available_start_positions)
        
        # Проверяем, что все позиции в батче мутабельны
        batch_positions = list(range(position, position + batch_size))
        if all(pos in mutable_positions for pos in batch_positions):
            # Выполняем скрещивание
            donor_part = donor_seq.upper()[position:position+batch_size]
            new_seq[position:position+batch_size] = list(donor_part)
    
    # Восстанавливаем фиксированные позиции из шаблона
    for pos, char in fixed_positions.items():
        if pos < len(new_seq):
            new_seq[pos] = char
    
    return "".join(new_seq)


def multiple_crossover_with_template(
    seq_list, 
    template="", 
    strict_length=False,
    max_mutations=10, 
    crossover_cycles=2, 
    max_crossover_batch=8, 
    material='dna'
):
    """
    Выполняет множественное скрещивание и мутацию с учетом шаблона.
    
    Args:
        seq_list: Список последовательностей
        template: Шаблон последовательности (A/T/G/C/U - фиксированные, x - переменные)
        strict_length: Если True, длина должна точно совпадать с шаблоном
        max_mutations: Максимальное количество мутаций
        crossover_cycles: Количество циклов скрещивания
        max_crossover_batch: Максимальный размер батча для скрещивания
        material: Тип материала ('dna' или 'rna')
        
    Returns:
        str: Результирующая последовательность, соответствующая шаблону
    """
    # Если шаблон не задан, используем обычную функцию
    if not template:
        return multiple_crossover(seq_list, max_mutations, crossover_cycles, max_crossover_batch, material)
    
    template_len = len(template)
    
    # Создаём комплементарную последовательность
    new_seq = compl(seq_list[0], material=material)
    
    # Применяем ограничения длины, если strict_length=True
    if strict_length:
        if len(new_seq) < template_len:
            # Дополняем случайными нуклеотидами
            new_seq += random_seq(template_len - len(new_seq), material)
        elif len(new_seq) > template_len:
            # Обрезаем до нужной длины
            new_seq = new_seq[:template_len]
    
    # Применяем фиксированные позиции из шаблона
    fixed_positions = _get_fixed_positions(template)
    new_seq_list = list(new_seq.upper())
    for pos, char in fixed_positions.items():
        if pos < len(new_seq_list):
            new_seq_list[pos] = char
    new_seq = "".join(new_seq_list)
    
    # Скрещивание с остальными последовательностями (с учетом шаблона)
    for seq in seq_list[1:]:
        compl_seq = compl(seq, material=material)
        new_seq = crossover_with_template(
            compl_seq, 
            new_seq, 
            template, 
            max_crossover_batch, 
            crossover_cycles
        )
    
    # Для strict_length=False шаблон применяется только к началу последовательности
    # Для strict_length=True шаблон применяется ко всей последовательности
    if strict_length:
        # Мутации с учетом шаблона (ко всей последовательности)
        mutable_positions = _get_mutable_positions(template)
        if mutable_positions:
            # Ограничиваем количество мутаций
            max_mutations = min(len(mutable_positions), max_mutations, len(new_seq))
            if max_mutations > 0:
                new_seq = mutate_x_letters_with_template(
                    new_seq, 
                    random.randint(1, max(1, max_mutations)), 
                    template, 
                    material
                )
    else:
        # strict_length=False: шаблон применяется только к первым len(template) символам
        # Если последовательность короче шаблона, дополняем
        if len(new_seq) < template_len:
            new_seq += random_seq(template_len - len(new_seq), material)
        
        # Применяем шаблон только к префиксу
        prefix = new_seq[:template_len]
        # Применяем фиксированные позиции к префиксу
        fixed_positions = _get_fixed_positions(template)
        prefix_list = list(prefix)
        for pos, char in fixed_positions.items():
            if pos < len(prefix_list):
                prefix_list[pos] = char
        prefix = "".join(prefix_list)
        
        # Мутации только в мутабельных позициях префикса
        mutable_positions = _get_mutable_positions(template)
        if mutable_positions:
            max_mutations = min(len(mutable_positions), max_mutations, len(prefix))
            if max_mutations > 0:
                prefix = mutate_x_letters_with_template(
                    prefix, 
                    random.randint(1, max(1, max_mutations)), 
                    template, 
                    material
                )
        
        # Объединяем префикс с остальной частью
        new_seq = prefix + new_seq[template_len:]
    
    return new_seq



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