import unittest

from lenrmc.koonin import LogLambda


class LogLambdaTest(unittest.TestCase):

    def test_p_p(self):
        rates = LogLambda(['p', 'p'], [1, 2, 5, 10]).rates()
        self.assertEqual(
            [-64.53426332253963,
             -48.17282760880237,
             -33.65482112413257,
             -26.337764631050703]
             , list(rates))

    def test_p_d(self):
        rates = LogLambda(['p', 'd'], [1, 2, 5, 10]).rates()
        self.assertEqual(
            [-54.98040222893245,
             -36.090262719865876,
             -19.32845937404553,
             -10.88053216470006]
             , list(rates))

    def test_p_t(self):
        rates = LogLambda(['p', 't'], [1, 2, 5, 10]).rates()
        self.assertEqual(
            [-57.762163335348745,
             -37.72974365009073,
             -19.954359876769544,
             -10.995589442761172]
             , list(rates))

    def test_d_d(self):
        rates = LogLambda(['d', 'd'], [1, 2, 5, 10]).rates()
        self.assertEqual(
            [-63.467401861202106,
             -40.33772257014951,
             -19.81404483116627,
             -9.470137792653718]
             , list(rates))

    def test_d_t(self):
        rates = LogLambda(['d', 't'], [1, 2, 5, 10]).rates()
        self.assertEqual(
            [-68.80088742375963,
             -43.472165963539396,
             -20.997210290845665,
             -9.669861697203626]
             , list(rates))
