import sqlite3


#search QUery: SELECT prime_extension.prime, prime_extension.k, l,  isogeny_eval.A, isogeny_eval.iso,  xADD, costs.xDBL, nr_add, costs.nr_mult, nr_div, nr_square FROM prime_extension, costs, isogeny_eval WHERE costs.id_prime = prime_extension.id_prime AND isogeny_eval.id_costs = costs.id_costs ORDER by prime;

class Database_iso():


    def __init__(self, p, k, l, costs, iso, A, file_name):
        self.p = int(p)
        self.k = int(k)
        self.costs = costs
        self.iso = iso
        self.l = int(l)
        self.A = int(A)
        self.conn = sqlite3.connect(file_name)


    @staticmethod
    def create_database(file_name):
        conn = sqlite3.connect(file_name)
        try:
            c = conn.cursor()
            table_1 = """ CREATE TABLE IF NOT EXISTS prime_extension (
                                        id_prime	INTEGER NOT NULL,
                                    	prime	INTEGER NOT NULL,
                                    	k	INTEGER NOT NULL,
                                    	PRIMARY KEY(id_prime AUTOINCREMENT)
                                    ); """
            c.execute(table_1)
        except Error as e:
            print(e)

        try:
            c = conn.cursor()
            table_2 = """ CREATE TABLE IF NOT EXISTS costs (
                                        id_costs	INTEGER NOT NULL,
                                    	id_prime	INTEGER NOT NULL,
                                    	frob	INTEGER,
                                    	xADD	INTEGER,
                                    	xDBL	INTEGER,
                                    	nr_mult	INTEGER,
                                    	nr_add	INTEGER,
                                    	nr_div	INTEGER,
                                    	nr_square	INTEGER,
                                    	l	INTEGER NOT NULL,
                                    	PRIMARY KEY(id_costs AUTOINCREMENT),
                                    	FOREIGN KEY(id_prime) REFERENCES prime_extension(id_prime)
                                    ); """
            c.execute(table_2)
        except Error as e:
            print(e)

        try:
            c = conn.cursor()
            table_3 = """ CREATE TABLE IF NOT EXISTS isogeny_eval (
                                        id_iso	INTEGER NOT NULL,
                                    	iso	BLOB,
                                    	id_costs	INTEGER NOT NULL,
                                    	A	INTEGER NOT NULL,
                                    	PRIMARY KEY(id_iso AUTOINCREMENT),
                                    	FOREIGN KEY(id_costs) REFERENCES costs(id_costs)
                                    ); """
            c.execute(table_3)
        except Error as e:
            print(e)





    def insert(self):
        print("starting insert...")
        c = self.conn.cursor()
        c.execute('INSERT INTO prime_extension (prime, k) VALUES (?, ?)', (self.p, self.k))
        prime_id = c.lastrowid
        self.conn.commit()


        self.__insert__costs__and__iso__(prime_id)

    def __insert__costs__and__iso__(self, prime_id):
        frob = 0;
        xADD = 0;
        xDBL = 0
        mult = 0
        add = 0
        div = 0
        square = 0
        if "frob" in self.costs:
            frob = self.costs["frob"]
        if "xADD" in self.costs:
            xADD = self.costs["xADD"]
        if "xDBL" in self.costs:
            xDBL = self.costs["xDBL"]
        if "mult" in self.costs:
            mult = self.costs["mult"]
        if "div" in self.costs:
            div = self.costs["div"]
        if "square" in self.costs:
            square = self.costs["square"]
        if "add" in self.costs:
            add = self.costs["add"]

        print(self.costs)
        c = self.conn.cursor()
        c.execute('INSERT INTO costs (id_prime, frob, xADD, xDBL, nr_mult, nr_add,  nr_div, nr_square, l) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)', (prime_id, frob, xADD, xDBL, mult, add,  div, square, self.l))
        c.execute('INSERT INTO isogeny_eval (iso, id_costs, A) VALUES (?, ?, ?)', (self.iso, c.lastrowid, self.A))
        self.conn.commit()