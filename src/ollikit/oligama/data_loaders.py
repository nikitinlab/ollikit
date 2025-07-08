import json
import numpy as np
import re
from pathlib import Path
from .exceptions import OligamaException, OligamaWarning
from .predictors.affinity_predictors import *
from .predictors.hairpin_predictors import *

class Dataloader():
	"""
    Базовый класс для загрузки и обработки данных.
    """
	def __init__(self, input_data, output_folder, seq_length=25):
		"""
        Инициализирует объект Dataloader.

        Args:
            input_data (str или dict): Путь к файлу JSON или словарь с данными.
            output_folder (str): Путь к папке для выходных данных.
            seq_length (int): Максимальная длина последовательности (по умолчанию 25).
        """     
		self.data = {}
		self.base_data_load(input_data, output_folder)
  
		self.seq_length = seq_length

		
		
	def add_init(self) :
     
		"""
        Инициализирует дополнительные атрибуты объекта.
        """     
        
		self.init_seqs = np.array(self.data["target_seqs"])
		self.target_seqs_orig = np.array(self.data["target_seqs"])
		self.target_seqs = self.prepare_clean_sequences()

		# target_names по умолчанию
		if "target_names" not in self.data or not self.data["target_names"]:
			self.target_names = [f"s{i+1}" for i in range(len(self.data["target_seqs"]))]
		else:
			self.target_names = self.data["target_names"]


		self.metric = self.data.get("metric", "fraction")

		self.celsius = float(self.data.get("celsius", 25))
		

		if self.celsius > 100 or self.celsius < 0:
			raise OligamaException("Temperature must be < 100 and > 0", self)
     
		self.aff_predictor_name = self.data.get("aff_predictor", "Nupack")
		self.hairpin_predictor_name = self.data.get("hairpin_predictor", "Seqfold")
 
		self.check_seqs_symbols()
  
		# Инициализация output_format
		self.output_format = self.data.get("output_format", "DNA").upper()
		#if self.output_format not in ["DNA", "RNA"]:
		#	raise OligamaException("output_format должен быть DNA или RNA", self)

		# Инициализация sequence_type
		self.sequence_types = self.data.get("sequence_types", ["DNA"] * len(self.data.get("target_seqs", [])))
		if not all(stype in ["DNA", "RNA"] for stype in self.sequence_types):
			raise OligamaException("All elements in sequence_type must be either 'DNA' or 'RNA'", self)

		#OligamaWarning(json.dumps(self.data, indent=4), self)
		#OligamaWarning(f"sequence_type: {set(self.sequence_types)}", self)
  
		# временное решение, т.к. не понимаю задачу. Предиктор получает одно значение meterial 
		if len(set(self.sequence_types)) > 1:
			raise OligamaException('All sequences must have the same type (DNA or RNA)', self)

		if len(self.sequence_types) > 0 and len(self.sequence_types) != len(self.data.get("target_seqs", [])):
			raise OligamaException("The length of sequence_type must match the number of sequences.", self)

		# Если все последовательности одного типа, устанавливаем sequences_type как этот тип
		unique_types = set(self.sequence_types)
		if len(unique_types) == 1:
			self.sequences_type_for_predictor = unique_types.pop()  # Устанавливаем один тип, если он единственный
		else:
			self.sequences_type_for_predictor = "DNA" # self.sequence_type  # Оставляем массив, если типы разные
  
		self.set_predictors()
  
		if "target_seqs" in self.data and self.aff_predictor_name == 'Oligama':
			self.check_seqs_len()
		elif "target_seqs" not in self.data:
			raise OligamaException("Input data must contain 'target_seqs'", self)

	def set_predictors(self):
     
		"""
        Устанавливает предикторы афинности и шпильки.
        """
		
		if self.aff_predictor_name == 'Oligama':
			self.aff_predictor = CNN_Affinity_Predictor(celsius = self.celsius)
		if self.hairpin_predictor_name == 'Oligama':
			self.hairpin_predictor = CNN_Hairpin_Predictor(celsius = self.celsius)
		if self.aff_predictor_name == 'Nupack':
			self.aff_predictor = Nupack_Affinity_Predictor(celsius = self.celsius, material=self.sequences_type_for_predictor.lower())
		if self.hairpin_predictor_name == 'Nupack':
			self.hairpin_predictor = Nupack_Hairpin_Predictor(celsius = self.celsius)
		if self.hairpin_predictor_name == 'Seqfold':
			self.hairpin_predictor = Seqfold_Hairpin_Predictor(celsius = self.celsius)

	def prepare_clean_sequences(self):
		"""
        Очищает последовательности от содержимого в скобках.

        Returns:
            np.array: Массив очищенных последовательностей.
        """
		cleaned_seqs = []
		for seq in self.target_seqs_orig:
			# Удаление содержимого в квадратных или круглых скобках
			cleaned_seq = re.sub(r"[\[\(].*?[\]\)]", "", seq)
			cleaned_seqs.append(cleaned_seq)
		return np.array(cleaned_seqs)

	def base_data_load(self, input_data, output_folder) :
		"""
        Загружает данные из файла или словаря.

        Args:
            input_data (str или dict): Путь к файлу JSON или словарь с данными.
            output_folder (str): Путь к папке для выходных данных.
        """
		self.output_folder = Path(output_folder)
		self.log_file = output_folder/Path("Log.txt")
		if isinstance(input_data, str):
            # Загружаем данные из файла
			with open(input_data) as f:
				self.data = json.load(f)
		elif isinstance(input_data, dict):
            # Используем переданные JSON-данные
			self.data = input_data
		else:
			raise OligamaException("input_data должен быть либо путем к файлу, либо JSON-объектом", self)
     
	def check_seqs_len(self):
		"""
        Проверяет длину последовательностей.
        """     
		len_arr = [len(seq) for seq in self.target_seqs]
		# проверка на одинаковость
		if len(set(len_arr)) > 1:
			raise OligamaException("Currently, all sequences must have the same length", self)

		for seq in self.target_seqs:
			if len(seq) > self.seq_length:
				raise OligamaException("Currently, the length of sequences must not exceed 25", self)


	def check_seqs_symbols(self):
		"""
        Проверяет символы в последовательностях.
        """     
		self.target_seqs = [seq.upper().replace(" ", "") for seq in self.target_seqs]
		for i in range(len(self.target_seqs)):
			self.target_seqs[i] = re.sub(r"([\(\[]).*?([\)\]])", "", self.target_seqs[i])
			if hasattr(self, 'aff_predictor_name') and self.aff_predictor_name == 'Oligama' and re.search(r'U', self.target_seqs[i]):
				mes = (
                    f"RNA sequence supplied in sequence {i}: {self.target_seqs[i]}. "
                    "U will be replaced with T, and sequence will be treated as DNA sequence. "
                    "Results will be less accurate. (Oligama only)"
                )
				OligamaWarning(mes, self)

				if re.search(r'T', self.target_seqs[i]): 
					raise OligamaException("Target sequence contains both T and U", self)
				
				self.target_seqs[i] = re.sub("U", "T", self.target_seqs[i])
    		
			if not self.target_seqs[i]:  
				raise OligamaException(f"Target sequence {self.target_names[i]} (index {i}) is empty", self)

			
			if re.fullmatch(r'[ATGCU]+', self.target_seqs[i]) is None:
				raise OligamaException("Target sequence contains forbidden symbols", self)

		self.target_seqs = np.array(self.target_seqs)

class Olig_Finder_Dataloader(Dataloader):
    
	def __init__(self, input_data, output_folder):
		"""
        Инициализирует объект Olig_Finder_Dataloader

        Args:
            input_data (str или dict): Путь к файлу JSON или словарь с данными.
            output_folder (str): Путь к папке для выходных данных.
		"""     
		super().__init__(input_data, output_folder)
		self.add_init()
		
		self.target_aff_low = np.array([float(i) for i in self.data["target_aff_low"]])
		self.target_aff_high = np.array([float(i) for i in self.data["target_aff_high"]])

		self.hairpin_en_thr = np.array(float(self.data["Hairpin_energy_thr"]))
		self.timeout = self.data.get("timeout", 300)  
		self.metric = self.data["metric"]
		self.num_oligos = int(self.data["num_oligos"])

		self.check_aff_thr()

	def check_aff_thr(self):
		"""
        Проверяет  правильность  заданных  порогов  афинности.
		"""     
		if not (len(self.target_aff_low) ==
				len(self.target_aff_high) == 
				len(self.target_seqs)):
				raise OligamaException("Target sequences array must have the same length as target threshold conditions", self)

		aff_arr_2d = np.concatenate((self.target_aff_low[:, None], self.target_aff_high[:, None]), axis = 1)
		new_aff_low = np.min(aff_arr_2d, axis = 1)
		new_aff_high = np.max(aff_arr_2d, axis = 1)

		if self.metric == 'fraction' and (new_aff_high.any() > 1 or new_aff_low.any() < 0):
			OligamaException("Affinity with fraction metric cannot be greater than 1 or less than 0", self)
		elif self.metric == 'Kd' and new_aff_low.any() < 0:
			OligamaException("Affinity with Kd metric cannot be less than 0", self)
		elif self.metric == 'gibbs' and  new_aff_high.any() > 0:
			OligamaException("Affinity with gibbs metric cannot be greater than 0", self)

		if not ((self.target_aff_low == new_aff_low).all() == True):
			OligamaWarning("Low and high affinity thresholds have been swaped for some sequences", self)

		self.target_aff_low, self.target_aff_high = new_aff_low, new_aff_high



class Concentration_Solver_Dataloader(Dataloader):
	def __init__(self, input_data, output_folder, conc_min=1e-15, conc_max=1e-2):
		"""
        Инициализирует  объект  Concentration_Solver_Dataloader  
        
        Args:
            input_data (str или dict): Путь к файлу JSON или словарь с данными.
            output_folder (str): Путь к папке для выходных данных.
            conc_min (float): Минимальная концентрация (по умолчанию 1e-15).
            conc_max (float): Максимальная концентрация (по умолчанию 1e-2).
		"""
		super().__init__(input_data, output_folder)
		self.add_init()
  
		self.conc_min = conc_min
		self.conc_max = conc_max
		self.total_conc = self.data["total_conc"]
		self.total_conc = [float(i) for i in self.total_conc]
		self.input_seq_name =  self.data["input_seq_name"]
		self.output_seq_name = self.data["output_seq_name"]
		self.input_seq_bounds = self.data["input_seq_bounds"]
		self.input_seq_bounds = [float(i) for i in self.input_seq_bounds]
		self.n_points = self.data["n_points"]

		if self.input_seq_bounds[0] > self.input_seq_bounds[1]:
			mes = "Lower input concentration bounds must be lower than the upper bound \n"
			raise OligamaException(mes, self)
		for conc in self.total_conc:
			if conc > self.conc_max or conc < self.conc_min:
				mes = "Initial concentration must be in range from 10 mM to 1 fM"
				mes = f"Initial concentration must be in range from {self.conc_min} to {self.conc_max}"
				raise OligamaException(mes, self)

class Complex_Solver_Dataloader(Dataloader):
    
	def __init__(self, input_data, output_folder):
		"""
        Инициализирует  объект  Complex_Solver_Dataloader

        Args:
            input_data (str или dict): Путь к файлу JSON или словарь с данными.
            output_folder (str): Путь к папке для выходных данных.
		"""
		super().__init__(input_data, output_folder)
		self.add_init()

		self.total_conc = self.data["total_conc"]
		self.total_conc = [float(i) for i in self.total_conc]


class Complement_Dataloader(Dataloader):  
    def __init__(self, input_data, output_folder):
        """
        Инициализирует  объект  Complement_Dataloader

        Args:
            input_data (str или dict): Путь к файлу JSON или словарь с данными.
            output_folder (str): Путь к папке для выходных данных.
        """        
        super().__init__(input_data, output_folder)

        # Проверяем, что target_seqs и target_names существуют
        if "target_seqs" not in self.data or "target_names" not in self.data:
            raise OligamaException("Input data must contain both 'target_seqs' and 'target_names'", self)
        
        self.target_seqs = np.array(self.data["target_seqs"])
        self.target_names = np.array(self.data["target_names"])
        self.sequence_types = self.data.get("sequence_types", ["DNA"] * len(self.data.get("target_seqs", [])))
        
        # Проверка соответствия длин списков последовательностей и имён
        if len(self.target_seqs) != len(self.target_names):
            raise OligamaException("Number of target sequences and names must match", self)
        
        # Выполняем проверки символов и на пустые последовательности
        self.check_seqs_symbols()


class MicroRNA_Dataloader(Dataloader):  
    
    def __init__(self, input_data, output_folder):
        """
        Инициализирует  объект  MicroRNA_Dataloader

        Args:
            input_data (str или dict): Путь к файлу JSON или словарь с данными.
            output_folder (str): Путь к папке для выходных данных.
        """        
        super().__init__(input_data, output_folder)

        self.aff_predictor_name = self.data.get("aff_predictor", "Oligama") or "Oligama"
        self.hairpin_predictor_name = self.data.get("hairpin_predictor", "Seqfold") or "Seqfold"
 
        self.check_seqs_symbols()
  
        self.set_predictors()
        
        # Проверяем, что target_seqs и target_names существуют
        if "target_seqs" not in self.data or "target_names" not in self.data:
            raise OligamaException("Input data must contain both 'target_seqs' and 'target_names'", self)
        
        self.target_seqs = np.array(self.data["target_seqs"])
        self.target_names = np.array(self.data["target_names"])
        
        # Проверка соответствия длин списков последовательностей и имён
        if len(self.target_seqs) != len(self.target_names):
            raise OligamaException("Number of target sequences and names must match", self)
        
        # Выполняем проверки символов и на пустые последовательности
        self.check_seqs_symbols()
