#!/usr/bin/python
from Tkinter import *
from ttk import *
import ionize

class Application(Frame):
    ions = []
    concentrations = []

    data = ionize.get_db()

    def add_ion(self):
        # ion_name = self.ion_entry.get()
        ion_name = self.db_tree.focus()
        #for multiple selection
        # ion_name = self.db_tree.selection()
        c = float(self.concentration_entry.get())
        ion = ionize.load_ion(ion_name)
        if ion and c:
            self.ions.append(ion)
            self.concentrations.append(c)
        # self.ion_list_display.insert(INSERT, str(ion.name)+
        # ' ' + str(c) + '\n')
        self.sol_tree.insert('', 'end', ion.name, text=str(ion.name),
                            values=c)

    def calc_solution(self):
        T = self.temp_entry.get()
        if T:
            sol = ionize.Solution(self.ions, self.concentrations, T=T)
        else:
            sol = ionize.Solution(self.ions, self.concentrations)
        self.solution_list_display.config(state=NORMAL)
        self.solution_list_display.delete(1.0, END)
        text = ''
        text += 'pH: {:.3g}\n'.format(sol.pH)
        text += 'Ionic strength: {:.3g} M\n'.format(sol.I)
        text += 'Conductivity: {:.3g} S/m\n'.format(sol.conductivity())
        text += 'KRF: {:.3g}\n'.format(sol.kohlrausch())
        text += 'Alberty: {:.3g}\n'.format(sol.alberty())
        text += 'Jovin: {:.3g}\n'.format(sol.jovin())

        text += '\nIon mobilities:\n---------------\n'
        for ion in sol.ions:
            text += '{}: {:.3g} m^2/V/s\n'.format(ion.name, ion.effective_mobility())

        self.solution_list_display.insert(INSERT, text)
        self.solution_list_display.config(state=DISABLED)

    def ion_search(self):
        searchstring = self.ion_search_entry.get()
        if searchstring:
            ions = ionize.search_ion(searchstring)
        else:
            ions = self.data.keys()
        ions=sorted(ions)
        self.db_tree_update(ions)


    def db_tree_update(self, items):
        for i in self.db_tree.get_children():
            self.db_tree.delete(i)
        for i in items:
            self.db_tree.insert('', 'end', i, text=i,
                                values=self.data[i][0:3])

    def database_tree(self, parent, data):
        self.db_tree = Treeview(parent)
        ysb = Scrollbar(self, orient='vertical', command=self.db_tree.yview)
        self.db_tree.configure(yscroll=ysb.set)
        self.db_tree['columns'] = ('z', 'pKa', 'mu')
        self.db_tree.heading('#0', text='Ion')
        self.db_tree.heading('z', text='Valance')
        self.db_tree.heading('pKa', text='pKa')
        self.db_tree.heading('mu', text='mobility')

        for item in sorted(data.keys()):
            self.db_tree.insert('', 'end', item, text=item,
                                values=data[item][0:3])
        self.db_tree.pack({"side": "left"})


    def solution_tree(self, parent):
        self.sol_tree = Treeview(parent)
        self.sol_tree['columns'] = ('conc')
        self.sol_tree.heading('#0', text='Ion')
        self.sol_tree.heading('conc', text='Concentration')
        self.sol_tree.pack({"side": "left"})

    def createWidgets(self):
        self.button_frame = Frame(self)

        self.searchstring = StringVar()
        self.searchstring.trace("w", lambda name, index, mode, sv=self.searchstring: self.ion_search())
        self.ion_search_entry = Entry(self.button_frame, textvariable=self.searchstring)
        self.ion_search_entry.pack({"side": "left"})

        vc = self.register(self.valid_conc, )
        self.concentration_entry = Entry(self.button_frame,validate='all', validatecommand=(vc, '%P'))
        self.concentration_entry.pack({"side": "left"})

        self.temp_entry = Entry(self.button_frame)
        self.temp_entry.pack({"side": "left"})

        self.add_ion_button = Button(self.button_frame)
        self.add_ion_button["text"] = "Add ion."
        self.add_ion_button["command"] = self.add_ion
        self.add_ion_button.pack({"side": "left"})

        self.calc_solution_button = Button(self.button_frame)
        self.calc_solution_button["text"] = "Calculate solution."
        self.calc_solution_button["command"] = self.calc_solution
        self.calc_solution_button.pack({"side": "left"})

        self.quit_button = Button(self.button_frame)
        self.quit_button["text"] = "Quit."
        self.quit_button["command"] =  self.quit
        self.quit_button.pack({"side": "left"})

        self.button_frame.pack()

        self.display_frame = Frame(self)
        self.database_tree(self.display_frame, ionize.get_db())
        self.solution_tree(self.display_frame)
        # self.ion_list_display = Text(self.display_frame)
        # self.ion_list_display.pack()
        self.display_frame.pack()
        self.solution_list_display = Text(self)
        self.solution_list_display.config(state=DISABLED)
        self.solution_list_display.pack()

    def valid_conc(self, string):
        try:
            n = float(string)
            if n>=0:
                return True
        except:
            if string in ['', '.']:
                return True
        return False

    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.pack()
        self.createWidgets()

root = Tk()
# img = PhotoImage(file='ionize_icon_v1.gif')
# root.tk.call('wm', 'iconphoto', root._w, img)
app = Application(master=root)
app.mainloop()
root.destroy()
