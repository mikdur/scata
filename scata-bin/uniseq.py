import pickle, sys

from collections.abc import MutableMapping, ValuesView


class UniseqIterator:
    def __init__(self, uniseq):
        self.uniseq = uniseq
        self.iterator = iter(self.uniseq.parts)
        
    def __next__(self):
        n=next(self.iterator)
        return self.uniseq[n]

    def __iter__(self):
        return self


class UniseqValuesView:
    def __init__(self, uniseq):
        self.uniseq = uniseq
        self.key_iter = iter(self.uniseq.keys())

    def __iter__(self):
        return self

    def __next__(self):
        k = next(self.key_iter)
        return self.uniseq[k]

class UniseqItemsView:
    def __init__(self, uniseq):
        self.uniseq = uniseq
        self.key_iter = iter(self.uniseq.keys())

    def __iter__(self):
        return self

    def __next__(self):
        k = next(self.key_iter)
        return (k, self.uniseq[k])
    

class UniseqDB(MutableMapping):
    def __init__(self, basename, mode="w", part_size=3000):
        self.basename = basename
        self.mode = mode
        self.part_size=part_size
        self.parts = {}
        self.partnum = 0
        self.cache = {}
        self.cache[self.partnum] = {}
        self.synced = True

        if mode == "w":
            try:
                self.parts = pickle.load(open(self.basename + "_parts.pick", "rb"))
            except IOError:
                pass
        else:
            self.parts = pickle.load(open(self.basename + "_parts.pick", "rb"))


    def keys(self):
        return self.parts.keys()

    

    def values(self):
        return UniseqValuesView(self)

    def __len__(self):
        return len(self.parts)

    def has_key(self, key):
        return key in self.parts

    def __contains__(self, key):
        return key in self.parts

    def get(self, key, default=None):
        if key in self.parts:
            return self[key]
        return default

    def __getitem__(self, key):
        value = self.cache[self.parts[key]].get(key)
        if value != None:
            return value
        if key not in self.parts:
            raise KeyError("UniseqDB: key %s not found" % (repr(key)))
        sys.stdout.write("UniseqDB: loading parition %s\n" % ( repr(self.parts[key])))
        sys.stdout.flush()
        self.cache[self.parts[key]] = pickle.load(open(self.basename + "_" + str(self.parts[key]) + ".pick", "rb"))
        value = self.cache[self.parts[key]][key]
        return value

    def __setitem__(self, key, value):
        if self.mode != "w":
            raise KeyError("Cannot update read-only UniseqDB")
        if key in self.parts:
            self.cache[self.parts[key]][key]=value
            return
        if len(self.cache[self.partnum]) > self.part_size:
            self.partnum += 1
            self.cache[self.partnum] = {}
        self.parts[key] = self.partnum
        self.cache[self.partnum][key] = value
        self.synced = False
        
    def __delitem__(self, key):
        if self.mode != "w":
            raise KeyError("Cannot update read-only UniseqDB")
        del self.cache[self.parts[key]][key]
        del self.parts[key]
        self.synced = False


    def close(self):
        if self.synced:
            return
        if self.mode == "w":
            self.sync()

    def __del__(self):
        self.close()

    def sync(self):
        if self.mode != "w":
            raise Exception("Will not sync a readonly db")
        sys.stdout.write("UniseqDB: syncing %d partitions\n" % (len(self.cache)))
        sys.stdout.flush()
        for part in self.cache:
            pickle.dump(self.cache[part], open(self.basename + "_" + str(part) + ".pick","wb"))
        pickle.dump(self.parts, open(self.basename + "_" + "parts.pick", "wb"))
        self.synced = True

    def items(self):
        return UniseqItemsView(self)
            
    def __iter__(self):
        return UniseqIterator(self)
    
