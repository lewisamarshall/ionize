#!/usr/bin/python
from Tkinter import *
from ttk import *
import ionize

class Application(Frame):
    ions = []
    concentrations = []

    data = ionize.get_db()

    def add_ion(self):
        ion_name = self.db_tree.focus()
        #for multiple selection
        # ion_name = self.db_tree.selection()
        c = float(self.concentration_entry.get())
        ion = ionize.load_ion(ion_name)

        self.ions.append(ion)
        self.concentrations.append(c)

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
        self.db_tree.heading('mu', text='Mobility')

        for item in sorted(data.keys()):
            self.db_tree.insert('', 'end', item, text=item,
                                values=data[item][0:3])
        self.db_tree.pack({"side": "left"})


    def solution_tree(self, parent):
        self.sol_tree = Treeview(parent)
        self.sol_tree['columns'] = ('conc')
        self.sol_tree.heading('#0', text='Ion')
        self.sol_tree.heading('conc', text='Concentration (M)')
        self.sol_tree.pack({"side": "left"})
        self.sol_tree.bind("<Double-1>", self.set_conc_popup)

    def set_concentration(self,item):
        c = float(self.c_box.get())
        self.sol_tree.set(item, column='conc', value=c)
        self.concentrations[self.ions.index(ionize.load_ion(item))] = c

    def set_conc_popup(self, event):
        item = self.sol_tree.identify('item',event.x,event.y)
        self.conc_pop = Toplevel()
        self.conc_pop.title('Set {} concentration.'.format(self.sol_tree.item(item,"text")))

        self.c_box = Entry(self.conc_pop)
        self.c_box.pack()

        button1 = Button(self.conc_pop, text="set", command=lambda: self.set_concentration(item))
        button1.pack()

        button2 = Button(self.conc_pop, text="Close", command=self.conc_pop.destroy)
        button2.pack()


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

        self.button_frame.pack()

        self.display_frame = Frame(self)
        self.database_tree(self.display_frame, ionize.get_db())
        self.solution_tree(self.display_frame)
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
app.master.title("ionize")
app.mainloop()
try:
    root.destroy()
except:
    pass
