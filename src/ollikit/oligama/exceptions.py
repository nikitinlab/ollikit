import re
import warnings

class OligamaException(Exception):
    """
    Пользовательское исключение для Oligama.

    Args:
        message (str): Сообщение об ошибке.
        obj (object): Объект, вызвавший исключение.
    """   
    def __init__(self, message, obj):
     
        message = re.sub(r'[\s]+', " ", message)
        super().__init__(message)
        with open(obj.log_file, "a") as f:
            f.write(message + "\n")

class OligamaWarning(Warning):
    """
    Пользовательское предупреждение для Oligama.

    Args:
        message (str): Сообщение предупреждения.
        obj (object): Объект, вызвавший предупреждение.
    """            
    def __init__(self, message, obj):

        message = re.sub(r'[\s]+', " ", message)
        warnings.warn(message)
        with open(obj.log_file, "a") as f:
            f.write(message + "\n")
