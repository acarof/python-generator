import iago
import os

class Analyser(iago.Analyser):
        def setup(self):
                self.path = os.getcwd()


if __name__ == '__main__':
        a = Analyser()
        a.run()
