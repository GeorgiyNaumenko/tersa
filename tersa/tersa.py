from .ersa import ERSA
from .dealer import Dealer
from .party import Party


class TERSA:

    def __int__(self, t: int, n: int):

        self.t = t
        self.n = n
        self.dealer = Dealer(t, n)
        self.parties = []
        self.started = False
        self.set_parties()
        self.ersa = None

    def set_parties(self):
        self.parties.clear()
        for i in range(self.n): self.parties.append(Party(i))

    def start_session(self):

        self.dealer.generate_params()

        for party in self.parties:
            party.set_partial(self.dealer)

        self.started = True

        self.ersa = ERSA(self.dealer.eta, self.dealer.gamma)

    def terminate_session(self):

        self.dealer = Dealer(self.t, self.n)
        self.set_parties()
        self.started = False
        self.ersa = None

    def search_party(self, idx):
        for i, party in enumerate(self.parties):
            if party.idx == idx:
                return i
        return -1

    def encrypt(self, message):
        if self.started:
            return self.ersa.encrypt(message)
        else:
            raise AttributeError("Session doesn't running currently. "
                                 "Start the new session for encrypting and decrypting.")

    def decrypt(self, message, parties_idx):
        partial_keys = []
        for party_idx in parties_idx:
            party = self.parties[party_idx]
            if party.idx == party_idx:
                partial_keys.append(party.partial_key)
            else:
                idx = self.search_party(party_idx)
                if idx > -1:
                    partial_keys.append(self.parties[idx].partial_key)
        common_key = self.dealer.reproduce_common_private_key(partial_keys)
        return self.ersa.decrypt(message, common_key)
