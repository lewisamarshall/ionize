#!/usr/bin/python
from Tkinter import *
from ttk import *
import ionize

class Application(Frame):
    ions = []
    concentrations = []

    def add_ion(self):
        ion_name = self.ion_entry.get()
        ion_name = self.db_tree.focus()
        #for multiple selection
        # ion_name = self.db_tree.selection()
        print ion_name
        c = float(self.concentration_entry.get())
        ion = ionize.load_ion(ion_name)
        if ion and c:
            self.ions.append(ion)
            self.concentrations.append(c)
        self.ion_list_display.insert(INSERT, str(ion.name)+
        ' ' + str(c) + '\n')

    def calc_solution(self):
        sol = ionize.Solution(self.ions, self.concentrations)
        self.solution_list_display.delete(1.0, END)
        text = ''
        text += 'pH: {:.3g}\n'.format(sol.pH)
        text += 'Ionic strength: {:.3g} M\n'.format(sol.I)
        text += 'Conductivity: {:.3g} S/m\n'.format(sol.conductivity())
        text += 'KRF: {:.3g}\n'.format(sol.kohlrausch())
        text += 'Alberty: {:.3g}\n'.format(sol.alberty())
        text += 'Jovin: {:.3g}\n'.format(sol.jovin())
        self.solution_list_display.insert(INSERT, text)



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
        self.db_tree.pack()


    def createWidgets(self):
        self.button_frame = Frame(self)

        self.ion_entry = Entry(self.button_frame)
        self.ion_entry.pack({"side": "left"})

        self.concentration_entry = Entry(self.button_frame)
        self.concentration_entry.pack({"side": "left"})

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
        # self.quit_button["fg"]   = "red"
        self.quit_button["command"] =  self.quit
        self.quit_button.pack({"side": "left"})

        self.button_frame.pack()

        self.display_frame = Frame(self)
        self.database_tree(self.display_frame, ionize.get_db())
        self.ion_list_display = Text(self.display_frame)
        self.ion_list_display.pack()
        self.solution_list_display = Text(self.display_frame)
        self.solution_list_display.pack()
        self.display_frame.pack()

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
