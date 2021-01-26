# -*- coding: utf-8 -*-

import tensorflow as tf


class OneHotPDGID(tf.keras.layers.Layer):
    def call(self, inputs):
        return tf.one_hot(tf.cast(tf.math.abs(inputs[:, :, -2]) > 11, tf.int32), 2)


class OneHotCharge(tf.keras.layers.Layer):
    def call(self, inputs):
        return tf.one_hot(tf.cast(inputs[:, :, -1] > 0, tf.int32), 2)


class OneHotYear(tf.keras.layers.Layer):
    def call(self, inputs):
        return tf.one_hot(tf.cast(inputs[:, 0] - 2016, tf.int32), 3)


class BlackBtag(tf.keras.layers.Layer):
    def call(self, inputs):
        return tf.keras.layers.Concatenate()([inputs[:, :, :-1], 0 * inputs[:, :, -1:]])


class BlackNBtag(tf.keras.layers.Layer):
    def call(self, inputs):
        return inputs[:, :-1]
