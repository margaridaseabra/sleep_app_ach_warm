"""Synthetic tests for pipelines/run_07_rem_theta.py (no animal data uploaded)."""
import importlib.util
import unittest
from pathlib import Path

import numpy as np

SOURCE = Path(__file__).resolve().parents[1] / "pipelines" / "run_07_rem_theta.py"
spec = importlib.util.spec_from_file_location("rem_theta", SOURCE)
rem = importlib.util.module_from_spec(spec)
spec.loader.exec_module(rem)


class RemThetaTests(unittest.TestCase):
    def test_sleep_scores_and_censoring(self):
        scores = np.r_[np.ones(30), np.full(50, 2.0), np.zeros(9), np.ones(40), np.full(30, 2.0)]
        converted = rem.reclassify_short_wake(scores)
        self.assertTrue(np.all(converted[80:89] == 15))
        bouts = rem.rem_bouts(converted)
        self.assertEqual(len(bouts), 2)
        self.assertEqual(bouts[0], (30, 80, 15, False))
        self.assertEqual(bouts[1], (129, 159, None, True))

    def test_theta_spectral_peak(self):
        fs = 610.3515625
        seconds = 60
        t = np.arange(int(fs * seconds)) / fs
        eeg = (np.sin(2 * np.pi * 8.0 * t) + 0.15 * np.sin(2 * np.pi * 2.0 * t)).astype(np.float32)
        out = rem.theta_features(eeg, fs, 0, seconds)
        self.assertGreater(len(out), 40)
        self.assertAlmostEqual(float(out["theta_peak_hz"].median()), 8.0, delta=0.4)
        self.assertGreater(float(out["theta_dominant"].mean()), 0.9)
        self.assertTrue(np.all(out["seconds_to_end_s"].to_numpy() >= 0))

    def test_non_theta_is_not_dominant(self):
        fs = 200
        t = np.arange(fs * 60) / fs
        eeg = np.sin(2 * np.pi * 2 * t)
        out = rem.theta_features(eeg, fs, 0, 60)
        self.assertLess(float(out["theta_dominant"].mean()), 0.1)

    def test_undefined_score_is_censored(self):
        scores = np.r_[np.ones(20), np.full(30, 2.0), [-1.0], np.ones(10)]
        bout = rem.rem_bouts(scores)[0]
        self.assertTrue(bout[3])

    def test_time_bins(self):
        self.assertEqual(rem.period_label(0), "ZT3_6")
        self.assertEqual(rem.period_label(10800), "ZT6_9")
        self.assertEqual(rem.period_label(21600), "ZT9_12")


if __name__ == "__main__":
    unittest.main()
