# app/server.py

import json
from http.server import BaseHTTPRequestHandler, HTTPServer
from urllib.parse import urlparse, parse_qs
from app.calculations import fetch_pdb_for_protein, start_simulation

class SimpleHTTPRequestHandler(BaseHTTPRequestHandler):
    def do_GET(self):
        parsed_path = urlparse(self.path)
        query_params = parse_qs(parsed_path.query)

        # Serve the HTML page (frontend)
        if parsed_path.path == "/":
            self.send_response(200)
            self.send_header("Content-type", "text/html")
            self.end_headers()
            try:
                with open("index.html", "r") as file:
                    self.wfile.write(file.read().encode())
            except FileNotFoundError:
                self.send_response(404)
                self.end_headers()
                self.wfile.write(b"Error: index.html not found.")

        elif parsed_path.path == "/fetch_pdb":
            protein_name = query_params.get("protein_name", [None])[0]
            ligand_name = query_params.get("ligand_name", [None])[0]
            
            if protein_name and ligand_name:
                # Fetch PDB for the protein
                pdb_content = fetch_pdb_for_protein(protein_name)
                if pdb_content:
                    response = {
                        "status": "success",
                        "message": f"PDB file for {protein_name} downloaded successfully."
                    }
                    self.send_response(200)
                    self.send_header("Content-type", "application/json")
                    self.end_headers()
                    self.wfile.write(json.dumps(response).encode())
                else:
                    self.send_response(500)
                    self.end_headers()
                    self.wfile.write(b"Failed to fetch PDB file.")
            else:
                self.send_response(400)
                self.end_headers()
                self.wfile.write(b"Error: Missing protein_name or ligand_name parameter.")
        
        elif parsed_path.path == "/run_simulation":
            protein_name = query_params.get("protein_name", [None])[0]
            
            if protein_name:
                # Start the simulation
                message = start_simulation(protein_name)
                response = {
                    "status": "success",
                    "message": message
                }
                self.send_response(200)
                self.send_header("Content-type", "application/json")
                self.end_headers()
                self.wfile.write(json.dumps(response).encode())
            else:
                self.send_response(400)
                self.end_headers()
                self.wfile.write(b"Error: Missing protein_name parameter.")
        
        else:
            self.send_response(404)
            self.end_headers()
            self.wfile.write(b"Page not found")

# Start the HTTP server
if __name__ == "__main__":
    server_address = ("", 8080)
    httpd = HTTPServer(server_address, SimpleHTTPRequestHandler)
    print("Server running on port 8080...")
    httpd.serve_forever()


def do_GET(self):
    parsed_path = urlparse(self.path)
    if parsed_path.path == '/run_simulation':
        params = parse_qs(parsed_path.query)
        protein_name = params.get('protein_name', [''])[0]
        ligand_name = params.get('ligand_name', [''])[0]
        
        if protein_name and ligand_name:
            message = start_simulation(protein_name, ligand_name)
            self.send_response(200)
            self.send_header('Content-type', 'text/plain')
            self.end_headers()
            self.wfile.write(message.encode())
        else:
            self.send_error(400, "Missing protein_name or ligand_name parameter")