from sys import exit, argv
from os import getcwd, listdir, environ, system, walk
from os.path import join, expanduser, exists, expandvars, isdir, splitext, basename
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QMainWindow, QFileDialog
from PyQt5.QtGui import QStandardItem, QStandardItemModel, QMovie
from PyQt5.uic import loadUi
from PIL import Image

class GIF_Generator(QMainWindow):
    def __init__(self):
        super().__init__()

        ui_location = join(getcwd(), "ui/gif_generator_layout.ui")
        loadUi(ui_location, self)

        self.btn_browse_for_folder.clicked.connect(self.browse_for_folder)
        self.val_opened_folder.returnPressed.connect(self.check_and_add_files)
        self.val_name_pattern.returnPressed.connect(self.filter_files)
        self.btn_sort_files.clicked.connect(self.sort_files)
        self.btn_create_gif.clicked.connect(self.create_gif)
    
    def browse_for_folder(self):
        default_directory = getcwd()
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        self.jpg_folder = QFileDialog.getExistingDirectory(self, "Browse Folder",
                                                default_directory, options=options)
        self.val_opened_folder.setText(self.jpg_folder)
        if self.jpg_folder:
            self.add_files_to_list()
            self.display_list(self.images)
    
    def check_and_add_files(self):
        try_path = self.val_opened_folder.text()
        if isdir(try_path):
            self.jpg_folder = try_path
            self.add_files_to_list()
            self.display_list(self.images)
        else:
            self.val_opened_folder.setText("Path does not exist")
    
    def add_files_to_list(self):        
        self.images = []
        for x in listdir(self.jpg_folder):
            if x.endswith(".png") or x.endswith(".jpg") or x.endswith(".jpeg"):
                self.images.append(x)
        self.filtered_images = self.images                      
    
    def display_list(self, chosen_list):
        self.list_files.clear()                
        for i in chosen_list:
            self.list_files.addItem(i)
        self.sort_files()             
    
    def filter_files(self):
        self.filtered_images = []
        text = self.val_name_pattern.text()
        for image in self.images:
            if text in image:
                self.filtered_images.append(image)        
        self.display_list(self.filtered_images)
    
    def sort_files(self):
        self.list_files.sortItems()
    
    def create_gif(self):
        output_gif_path = join(self.jpg_folder,"output.gif")
        self.filtered_images.sort()
        image_paths = [join(self.jpg_folder,x) for x in self.filtered_images]
        images = [Image.open(image_path) for image_path in image_paths]
        images[0].save(output_gif_path,save_all=True,append_images=images[1:],duration=250,loop=0)
        self.movie = QMovie(output_gif_path)
        self.image_output.setMovie(self.movie)
        self.movie.start()


        
    

    

    




        

if __name__ == "__main__":
    app = QApplication(argv)
    window = GIF_Generator()
    window.show()
    exit(app.exec_())



