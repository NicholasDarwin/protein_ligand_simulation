import http.server
import socketserver
import json
from urllib.parse import urlparse, parse_qs
from app.calculations import fetch_pdb_for_protein, start_simulation

PORT = 8080

class MyHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    def do_GET(self):
        parsed_path = urlparse(self.path)
        path = parsed_path.path
        query_params = parse_qs(parsed_path.query)

        if path == "/":
            try:
                with open("index.html", 'r') as f:
                    html_content = f.read()
                self.send_response(200)
                self.send_header('Content-type', 'text/html')
                self.end_headers()
                self.wfile.write(html_content.encode())
            except FileNotFoundError:
                self.send_error(404, "File Not Found")

        elif path == "/fetch_pdb":
            protein_name = query_params.get('protein_name', [''])[0]
            ligand_name = query_params.get('ligand_name', [''])[0]
            
            if protein_name and ligand_name:
                try:
                    result = fetch_pdb_for_protein(protein_name)
                    response = {
                        'success': True,
                        'message': f"PDB files fetched successfully for protein: {protein_name} and ligand: {ligand_name}"
                    }
                except Exception as e:
                    response = {
                        'success': False,
                        'message': f"Error fetching PDB: {str(e)}"
                    }
            else:
                response = {
                    'success': False,
                    'message': "Both protein and ligand names are required"
                }
            
            self.send_response(200)
            self.send_header('Content-type', 'application/json')
            self.end_headers()
            self.wfile.write(json.dumps(response).encode())

        elif path == "/run_simulation":
            protein_name = query_params.get('protein_name', [''])[0]
            ligand_name = query_params.get('ligand_name', [''])[0]

            if protein_name and ligand_name:
                try:
                    # Pass both protein_name and simulation_name to start_simulation
                    result = start_simulation(protein_name, ligand_name)
                    response = {
                        'success': True,
                        'message': result
                    }
                except Exception as e:
                    response = {
                        'success': False,
                        'message': f"Error running simulation: {str(e)}"
                    }

            self.send_response(200)
            self.send_header('Content-type', 'application/json')
            self.end_headers()
            self.wfile.write(json.dumps(response).encode())

        else:
            self.send_error(404, "Endpoint not found")

if __name__ == "__main__":
    try:
        with socketserver.TCPServer(("", PORT), MyHTTPRequestHandler) as httpd:
            print(f"Server started at http://localhost:{PORT}")
            httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nServer stopped by user")
    except Exception as e:
        print(f"Error starting server: {e}")