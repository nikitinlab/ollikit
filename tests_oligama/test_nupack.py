import logging
import numpy as np
from ollikit import Nupack_Affinity_Predictor

# Настройка логирования
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('nupack_test.log'),
        logging.StreamHandler()
    ]
)

def test_complementary_sequences():
    # Комплементарные последовательности
    seq1 = "ATCGATCGATCG"
    seq2 = "TAGCTAGCTAGC"  # Комплементарная к seq1
    
    # Не комплементарная последовательность
    seq3 = "GGGGGGGGGGGG"
    
    predictor = Nupack_Affinity_Predictor(celsius=25)
    
    # Тест комплементарных последовательностей
    logging.info("\nТест комплементарных последовательностей:")
    try:
        affinity = predictor.simple_aff(seq1, seq2, units='fraction')
        logging.info(f"Аффинность между комплементарными последовательностями: {affinity:.4f}")
    except Exception as e:
        logging.error(f"Ошибка при расчете аффинности: {e}")
    
    # Тест не комплементарных последовательностей
    logging.info("\nТест не комплементарных последовательностей:")
    try:
        affinity = predictor.simple_aff(seq1, seq3, units='fraction')
        logging.info(f"Аффинность между не комплементарными последовательностями: {affinity:.4f}")
    except Exception as e:
        logging.error(f"Ошибка при расчете аффинности: {e}")

if __name__ == "__main__":
    test_complementary_sequences() 