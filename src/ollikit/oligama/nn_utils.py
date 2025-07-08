import keras
import tensorflow as tf
import tensorflow.keras.saving
import tensorflow.keras.backend as K

import pandas as pd

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
