from .dealer import Dealer


class Party:

    def __init__(self, idx: int):
        self.partial_key = None
        self.idx = idx

    def set_partial(self, dealer: Dealer):
        self.partial_key = dealer.partial_secrets[self.idx]

    def get_partial(self):
        return self.partial_key
