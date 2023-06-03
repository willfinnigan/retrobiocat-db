""" Module containing the LSTM-based model for calculation route distances """
from typing import List, Tuple, Dict, Any

import torch
import pytorch_lightning as lightning
from treelstm import TreeLSTM as TreeLSTMBase
from torchmetrics import MeanAbsoluteError, R2Score

from retrobiocat_web.retro.pathway_search.route_evaluation.clustering.route_distances import defaults

StrDict = Dict[str, Any]

class _TreeLstmWithPreCompression(torch.nn.Module):
    def __init__(self, fp_size: int, lstm_size: int, dropout_prob: float) -> None:
        super().__init__()
        self._compression = torch.nn.Sequential(
            torch.nn.Linear(fp_size, lstm_size),
            torch.nn.ReLU(),
            torch.nn.Dropout(p=dropout_prob),
            torch.nn.Linear(lstm_size, lstm_size),
            torch.nn.ReLU(),
        )
        self._tree_lstm = TreeLSTMBase(lstm_size, lstm_size)

    def forward(self, tree_batch: StrDict) -> torch.Tensor:
        """
        Forward pass
        :param tree_batch: collated trees from the `route_distances.utils.collate_trees` function.
        :return: the LSTM representation of the first nodes
        """
        features = self._compression(tree_batch["features"])
        lstm_output, _ = self._tree_lstm(
            features,
            tree_batch["node_order"],
            tree_batch["adjacency_list"],
            tree_batch["edge_order"],
        )
        # Only save value of top-node
        lstm_output = torch.stack(
            [t[0, :] for t in torch.split(lstm_output, tree_batch["tree_sizes"], dim=0)]
        )
        return lstm_output


class RouteDistanceModel(lightning.LightningModule):
    """
    Model for computing the distances between two synthesis routes
    :param fp_size: the length of the fingerprint vector
    :param lstm_size: the size o the LSTM cell
    :param dropout_prob: the dropout probability
    :param learning_rate: the initial learning rate of the optimizer
    :param weight_decay: weight decay factor of the optimizer
    """

    def __init__(
        self,
        fp_size: int = defaults.FP_SIZE,
        lstm_size: int = defaults.LSTM_SIZE,
        dropout_prob: float = defaults.DROPOUT_PROB,
        learning_rate: float = defaults.LEARNING_RATE,
        weight_decay: float = defaults.WEIGHT_DECAY,
    ) -> None:
        super().__init__()
        self.save_hyperparameters()
        self._tree_lstm = _TreeLstmWithPreCompression(fp_size, lstm_size, dropout_prob)
        self._pdist = torch.nn.PairwiseDistance(p=2)
        self._loss_func = torch.nn.MSELoss()
        self._mae = MeanAbsoluteError()
        self._r2 = R2Score()
        self._lr = learning_rate
        self._weight_decay = weight_decay

    def forward(self, tree_data: StrDict) -> torch.Tensor:
        """
        Calculate the pairwise distances between the input trees
        :param tree_data: collated trees from the `route_distances.utils.collate_trees` function.
        :return: the distances in condensed form
        """
        lstm_enc = self._tree_lstm(tree_data)
        return torch.pdist(lstm_enc)
