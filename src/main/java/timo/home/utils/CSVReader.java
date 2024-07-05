package timo.home.utils;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

public class CSVReader{
	public ArrayList<ArrayList<String>> data;
	
	public CSVReader(String fileName){
		this(fileName,",");
	}

	public CSVReader(String fileName,String separator){
		data = new ArrayList<ArrayList<String>>();
		try {
			BufferedReader br = new BufferedReader(new FileReader(fileName));
			String line;
			while ((line = br.readLine()) != null) {
				String[] split = line.split(separator);
				data.add(new ArrayList<String>());
				for (int i =0;i<split.length;++i){
					data.get(data.size()-1).add(split[i]);
				}
			}
		} catch (Exception e) {System.out.println("Couldn't read the file");}
	}

	public void printData(){
		System.out.println("Print data");	
		for (int i =0;i<data.size();++i){
			for (int j =0;j<data.get(i).size();++j){
				System.out.print(data.get(i).get(j)+"\t");
			}	
			System.out.println("");	
		}
		System.out.println("");	
	}
}
