#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  3 12:14:45 2025

@author: luismoncayo
"""
class UploadDataLineBalacing:
    
    def __init__(self, data_path):
        self.path_to_instance = data_path
        self.tasks = []
        self.tasks_times =[]
        self.precedences = []
        self.inst_name = "name"
        self.input_cycle_time = 0
    
    def upload_data(self):
        print(self.path_to_instance)
        parts = self.path_to_instance.split("instance", 1)
        model_name = "None"
        if len(parts) > 1:
            result = parts[1].strip()  # Get the part after the keyword
            after_dot = result.split('.', 1)[0]
            #print("instance"+after_dot)  # Output: "jumps over the lazy dog"
            model_name = "instance"+after_dot
            self.inst_name = model_name
        else:
            print("instance word is not found")
        
        with open(self.path_to_instance, 'r') as f:
            lines = [line.rstrip('\n') for line in f]
            #print(lines)
        f.close()
        
        number_tasks = int(lines[lines.index('<number of tasks>')+1])
        cycle_time =  int(lines[lines.index('<cycle time>')+1])
        self.input_cycle_time = cycle_time
        
        index_task_time = lines.index('<task times>')
        tasks_id = []
        tasks_times = []
        for i in range(index_task_time+1, index_task_time+number_tasks+1):
            number = int(lines[i].split(" ")[0])
            time = int(lines[i].split(" ")[1])
            tasks_id.append(number)
            tasks_times.append(time)
        
        index_precedence = lines.index('<precedence relations>')
        precedence = []
        k = index_precedence+1
        while len(lines[k]) != 0:
            from_t = int(lines[k].split(",")[0])
            to_t = int(lines[k].split(",")[1])
            precedence.append((from_t, to_t))
            k+=1
        #print(tasks_num)
        self.tasks = tasks_id
        #print(tasks_times)
        self.tasks_times = tasks_times
        #print(precedence)
        self.precedences = precedence
        
        if not all(self.tasks.count(x) == 1 for x in self.tasks) :
            print("The task ID in the file is not correctly enumerated")
        sorted_lst = sorted(self.tasks)
        test_tasks = list(range(1,number_tasks+1))
        if test_tasks==sorted_lst:
            print()
            print(f"There are {len(sorted_lst)} consecutively numbered tasks.")
        else:
            print("The tasks are not enumerated consecutively.")
        print()
        if (number_tasks == len(self.tasks)) and (len(self.tasks) == len(self.tasks_times)):
            print("The vectors of tasks and the tasks times are of the same length.")
            print()
        else:
            print("The vectors of tasks and the tasks times are not of the same length,")
            print("or/and the number of task quoted in the read file is different.")
            print()
        for i in self.precedences:
            if (i[0] not in self.tasks) or (i[1] not in self.tasks) or (i[0] > i[1]):
                print("A precedence value is not in the tasks vector or is enumerated wrongly.")
        print("The precedence relationships (i,j) are labelled as i < j.") 
        print()
        
        return self.tasks,self.tasks_times,self.precedences
        
    def get_tasks(self):
        return self.tasks
    
    def get_tasks_times(self):
        return self.tasks_times
    
    def get_precedences(self):
        return self.precedences
    
    def get_instance_name(self):
        return self.inst_name
    
    def get_input_cycle_time(self):
        return self.input_cycle_time
    
    

