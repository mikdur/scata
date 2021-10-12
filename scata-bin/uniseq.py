import pickle, sys

import UserDict

class UniseqDB(UserDict.DictMixin):
    def __init__(self, basename, mode="w", part_size=3000):
        self.basename = basename
        self.mode = mode
        self.part_size=part_size
        self.parts = {}
        self.partnum = 0
        self.cache = {}
        self.cache[self.partnum] = {}
        self.key_list = False
        self.synced = True

        if mode == "w":
            try:
                self.parts = pickle.load(open(self.basename + "_parts.pick"))
            except IOError:
                pass
        else:
            self.parts = pickle.load(open(self.basename + "_parts.pick"))


    def keys(self):
        if not self.key_list:
            self.key_list = list(self.parts.keys())
        return self.key_list

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
        try:
            value = self.cache[self.parts[key]][key]
        except KeyError:
            if key not in self.parts:
                raise KeyError("UniseqDB: key %s not found" % (repr(key)))
            sys.stdout.write("UniseqDB: loading parition %s\n" % ( repr(self.parts[key])))
            sys.stdout.flush()
            self.cache[self.parts[key]] = pickle.load(open(self.basename + "_" + str(self.parts[key]) + ".pick"))
            value = self.cache[self.parts[key]][key]
        return value

    def __setitem__(self, key, value):
        if key in self.parts:
            self.cache[self.parts[key]][key]=value
            return
        if len(self.cache[self.partnum]) > self.part_size:
            self.partnum += 1
            self.cache[self.partnum] = {}
        self.parts[key] = self.partnum
        self.cache[self.partnum][key] = value
        if not self.key_list:
            self.key_list = [ ]
        self.key_list.append(key)
        self.synced = False
        
    def __delitem__(self, key):
        del self.cache[self.parts[key]][key]
        del self.parts[key]
        del self.key_list[self.key_list.index(key)]
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
            pickle.dump(self.cache[part], open(self.basename + "_" + str(part) + ".pick","wct"))
        pickle.dump(self.parts, open(self.basename + "_" + "parts.pick", "wct"))
        self.synced = True
                   
    
